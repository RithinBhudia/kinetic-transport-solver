#include "kinetic.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define AT(s, i, j) ((s)->f[(size_t)(i) * (s)->cfg.nv + (j)])
#define STAGE(s, i, j) ((s)->stage[(size_t)(i) * (s)->cfg.nv + (j)])

/*
  The reference BLAS/LAPACK ABI exposes Fortran routines with a trailing
  underscore.  Calling DGTSV directly keeps this project usable with Apple's
  Accelerate, Netlib LAPACK, and vendor LAPACKs without requiring LAPACKE.
 */
extern void dgtsv_(const int *n, const int *nrhs, double *dl, double *d,
                   double *du, double *b, const int *ldb, int *info);

static double gaussian(double value, double center, double width) {
    double z = (value - center) / width;
    return exp(-0.5 * z * z);
}

int simulation_init(Simulation *s, const SimulationConfig *cfg) {
    size_t count;
    memset(s, 0, sizeof *s); s->cfg = *cfg;
    s->dx = (cfg->x_max - cfg->x_min) / (cfg->nx - 1);
    s->dv = (cfg->v_max - cfg->v_min) / (cfg->nv - 1);
    count = (size_t)cfg->nx * cfg->nv;
    s->x = malloc((size_t)cfg->nx * sizeof *s->x); s->v = malloc((size_t)cfg->nv * sizeof *s->v);
    s->force = malloc((size_t)cfg->nx * sizeof *s->force); s->f = calloc(count, sizeof *s->f); s->stage = calloc(count, sizeof *s->stage);
    if (!s->x || !s->v || !s->force || !s->f || !s->stage) { fprintf(stderr, "allocation failed\n"); simulation_destroy(s); return 1; }
    for (int i = 0; i < cfg->nx; ++i) {
        s->x[i] = cfg->x_min + i * s->dx;
        s->force[i] = cfg->force_offset + cfg->force_slope * (s->x[i] - 0.5 * (cfg->x_min + cfg->x_max));
    }
    for (int j = 0; j < cfg->nv; ++j) s->v[j] = cfg->v_min + j * s->dv;
    for (int i = 0; i < cfg->nx; ++i) for (int j = 0; j < cfg->nv; ++j)
        AT(s, i, j) = cfg->initial_amplitude * gaussian(s->x[i], cfg->initial_x_center, cfg->initial_x_width) * gaussian(s->v[j], cfg->initial_v_center, cfg->initial_v_width);
    return 0;
}

void simulation_destroy(Simulation *s) { free(s->x); free(s->v); free(s->force); free(s->f); free(s->stage); memset(s, 0, sizeof *s); }

static void choose_timestep(Simulation *s) {
    double v_max = fmax(fabs(s->cfg.v_min), fabs(s->cfg.v_max)), force_max = 0.0;
    for (int i = 0; i < s->cfg.nx; ++i) force_max = fmax(force_max, fabs(s->force[i]));
    double transport_limit = INFINITY;
    if (v_max > 0.0) transport_limit = fmin(transport_limit, s->dx / v_max);
    if (force_max > 0.0) transport_limit = fmin(transport_limit, s->dv / force_max);
    s->steps = s->cfg.min_steps > 1 ? s->cfg.min_steps : 1;
    if (isfinite(transport_limit)) { int cfl_steps = (int)ceil(s->cfg.final_time / (s->cfg.cfl * transport_limit)); if (cfl_steps > s->steps) s->steps = cfl_steps; }
    s->dt = s->cfg.final_time / s->steps;
}

static void apply_boundaries(Simulation *s) {
    int nx = s->cfg.nx, nv = s->cfg.nv;
    for (int j = 0; j < nv; ++j) {
        if (s->v[j] > 0.0) AT(s, 0, j) = s->cfg.left_inflow_amplitude * gaussian(s->v[j], s->cfg.left_inflow_center, s->cfg.left_inflow_width);
        else AT(s, 0, j) = AT(s, 1, j); /* zero-gradient outflow */
        if (s->v[j] < 0.0) AT(s, nx - 1, j) = s->cfg.right_inflow_amplitude * gaussian(s->v[j], s->cfg.right_inflow_center, s->cfg.right_inflow_width);
        else AT(s, nx - 1, j) = AT(s, nx - 2, j);
    }
    // Velocity Dirichlet data takes precedence where the boundaries meet. 
    for (int i = 0; i < nx; ++i) AT(s, i, 0) = AT(s, i, nv - 1) = 0.0;
}

static void advect_upwind(Simulation *s) {
    int nx = s->cfg.nx, nv = s->cfg.nv;
    for (int i = 0; i < nx; ++i) STAGE(s, i, 0) = STAGE(s, i, nv - 1) = 0.0;
    for (int j = 0; j < nv; ++j) { STAGE(s, 0, j) = AT(s, 0, j); STAGE(s, nx - 1, j) = AT(s, nx - 1, j); }
    for (int i = 1; i < nx - 1; ++i) for (int j = 1; j < nv - 1; ++j) {
        double dfdx = s->v[j] >= 0.0 ? (AT(s,i,j)-AT(s,i-1,j))/s->dx : (AT(s,i+1,j)-AT(s,i,j))/s->dx;
        double dfdv = s->force[i] >= 0.0 ? (AT(s,i,j)-AT(s,i,j-1))/s->dv : (AT(s,i,j+1)-AT(s,i,j))/s->dv;
        STAGE(s,i,j) = AT(s,i,j) - s->dt * (s->v[j] * dfdx + s->force[i] * dfdv);
    }
}

static int diffuse_implicit(Simulation *s) {
    int n = s->cfg.nv - 2;
    if (s->cfg.diffusion == 0.0) { memcpy(s->f, s->stage, (size_t)s->cfg.nx*s->cfg.nv*sizeof *s->f); return 0; }
    double r = s->dt * s->cfg.diffusion / (s->dv * s->dv);
    double *dl = malloc((size_t)(n-1)*sizeof *dl), *d = malloc((size_t)n*sizeof *d), *du = malloc((size_t)(n-1)*sizeof *du), *rhs = malloc((size_t)n*sizeof *rhs);
    if (!dl || !d || !du || !rhs) { free(dl); free(d); free(du); free(rhs); return 1; }
    for (int i = 0; i < s->cfg.nx; ++i) {
        for (int k = 0; k < n-1; ++k) dl[k] = du[k] = -r;
        for (int k = 0; k < n; ++k) { d[k] = 1.0 + 2.0*r; rhs[k] = STAGE(s,i,k+1); }
        int one = 1, info = 0;
        dgtsv_(&n, &one, dl, d, du, rhs, &n, &info);
        if (info != 0) { free(dl); free(d); free(du); free(rhs); return 1; }
        AT(s,i,0) = AT(s,i,s->cfg.nv-1) = 0.0;
        for (int k = 0; k < n; ++k) AT(s,i,k+1) = rhs[k];
    }
    free(dl); free(d); free(du); free(rhs); return 0;
}

int simulation_run(Simulation *s) {
    choose_timestep(s);
    printf("Running %d steps with dt = %.6e (CFL target %.3f)\n", s->steps, s->dt, s->cfg.cfl);
    for (int step = 0; step < s->steps; ++step) {
        apply_boundaries(s);
        advect_upwind(s);
        if (diffuse_implicit(s)) {
            fprintf(stderr, "LAPACK tridiagonal solve failed\n");
            return 1;
        }
        /* The implicit row solve changes x-edge rows, so restore the
         * characteristic data before the next step and before output. */
        apply_boundaries(s);
    }
    return 0;
}

int simulation_write_csv(const Simulation *s, const char *path) {
    FILE *out = fopen(path, "w"); if (!out) { perror(path); return 1; }
    fprintf(out, "x,v,f\n");
    for (int i = 0; i < s->cfg.nx; ++i) for (int j = 0; j < s->cfg.nv; ++j) fprintf(out, "%.16e,%.16e,%.16e\n", s->x[i], s->v[j], AT(s,i,j));
    return fclose(out) != 0;
}
