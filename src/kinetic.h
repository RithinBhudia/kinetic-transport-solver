#ifndef KINETIC_H
#define KINETIC_H

#include <stddef.h>

typedef struct {
    int nx, nv;
    double x_min, x_max, v_min, v_max;
    double final_time, cfl, diffusion;
    int min_steps;
    double force_offset, force_slope;
    double initial_amplitude, initial_x_center, initial_x_width;
    double initial_v_center, initial_v_width;
    double left_inflow_amplitude, left_inflow_center, left_inflow_width;
    double right_inflow_amplitude, right_inflow_center, right_inflow_width;
    char output_file[256];
} SimulationConfig;

typedef struct {
    SimulationConfig cfg;
    double dx, dv, dt;
    int steps;
    double *x, *v, *force;
    double *f, *stage;
} Simulation;

int config_read(const char *path, SimulationConfig *cfg);
int simulation_init(Simulation *sim, const SimulationConfig *cfg);
void simulation_destroy(Simulation *sim);
int simulation_run(Simulation *sim);
int simulation_write_csv(const Simulation *sim, const char *path);

#endif
