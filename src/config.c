#include "kinetic.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void defaults(SimulationConfig *c) {
    *c = (SimulationConfig){
        .nx = 101, .nv = 121, .x_min = 0.0, .x_max = 1.0,
        .v_min = -4.0, .v_max = 4.0, .final_time = 0.5,
        .cfl = 0.45, .diffusion = 0.05, .min_steps = 1,
        .force_offset = 0.0, .force_slope = 0.0,
        .initial_amplitude = 1.0, .initial_x_center = 0.5,
        .initial_x_width = 0.12, .initial_v_center = 0.0,
        .initial_v_width = 0.8,
        .left_inflow_amplitude = 0.0, .left_inflow_center = 0.0,
        .left_inflow_width = 1.0, .right_inflow_amplitude = 0.0,
        .right_inflow_center = 0.0, .right_inflow_width = 1.0,
    };
    strcpy(c->output_file, "output/solution.csv");
}

static int set_int(const char *key, const char *value, SimulationConfig *c) {
    if (!strcmp(key, "nx")) c->nx = atoi(value);
    else if (!strcmp(key, "nv")) c->nv = atoi(value);
    else if (!strcmp(key, "min_steps")) c->min_steps = atoi(value);
    else return 0;
    return 1;
}

static int set_double(const char *key, const char *value, SimulationConfig *c) {
    double v = strtod(value, NULL);
#define SET(name) if (!strcmp(key, #name)) c->name = v
    SET(x_min); else SET(x_max); else SET(v_min); else SET(v_max);
    else SET(final_time); else SET(cfl); else SET(diffusion);
    else SET(force_offset); else SET(force_slope);
    else SET(initial_amplitude); else SET(initial_x_center); else SET(initial_x_width);
    else SET(initial_v_center); else SET(initial_v_width);
    else SET(left_inflow_amplitude); else SET(left_inflow_center); else SET(left_inflow_width);
    else SET(right_inflow_amplitude); else SET(right_inflow_center); else SET(right_inflow_width);
    else return 0;
#undef SET
    return 1;
}

int config_read(const char *path, SimulationConfig *cfg) {
    char line[512], *key, *value, *comment;
    FILE *in;
    defaults(cfg);
    if (!(in = fopen(path, "r"))) { perror(path); return 1; }
    for (int line_no = 1; fgets(line, sizeof line, in); ++line_no) {
        comment = strchr(line, '#'); if (comment) *comment = '\0';
        key = line; while (isspace((unsigned char)*key)) ++key;
        if (*key == '\0') continue;
        value = strchr(key, '=');
        if (!value) { fprintf(stderr, "%s:%d: expected key = value\n", path, line_no); fclose(in); return 1; }
        *value++ = '\0';
        while (isspace((unsigned char)*value)) ++value;
        char *end = key + strlen(key); while (end > key && isspace((unsigned char)end[-1])) *--end = '\0';
        end = value + strlen(value); while (end > value && isspace((unsigned char)end[-1])) *--end = '\0';
        if (!strcmp(key, "output_file")) snprintf(cfg->output_file, sizeof cfg->output_file, "%s", value);
        else if (!set_int(key, value, cfg) && !set_double(key, value, cfg)) {
            fprintf(stderr, "%s:%d: unknown setting '%s'\n", path, line_no, key); fclose(in); return 1;
        }
    }
    fclose(in);
    if (cfg->nx < 3 || cfg->nv < 3 || cfg->x_max <= cfg->x_min || cfg->v_max <= cfg->v_min ||
        cfg->final_time < 0.0 || cfg->cfl <= 0.0 || cfg->diffusion < 0.0 || cfg->initial_x_width <= 0.0 ||
        cfg->initial_v_width <= 0.0 || cfg->left_inflow_width <= 0.0 || cfg->right_inflow_width <= 0.0) {
        fprintf(stderr, "invalid grid, physical parameter, or Gaussian width\n"); return 1;
    }
    return 0;
}
