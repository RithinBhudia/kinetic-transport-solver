#include "kinetic.h"
#include <stdio.h>

int main(int argc, char **argv) {
    SimulationConfig cfg; Simulation sim;
    const char *config_path = argc > 1 ? argv[1] : "examples/beam.cfg";
    if (argc > 2) { fprintf(stderr, "Usage: %s [config-file]\n", argv[0]); return 2; }
    if (config_read(config_path, &cfg) || simulation_init(&sim, &cfg)) return 1;
    int status = simulation_run(&sim);
    if (!status) status = simulation_write_csv(&sim, cfg.output_file);
    if (!status) printf("Wrote %s\n", cfg.output_file);
    simulation_destroy(&sim); return status;
}
