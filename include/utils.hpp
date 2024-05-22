#pragma once

#include <gsl/gsl_integration.h>

struct gsl_workspace {
    gsl_integration_workspace *ws;
    gsl_integration_qaws_table *table;
};

void free_ws(gsl_workspace ws);
/* gsl_workspace setup_ws(double alpha, double beta); */
gsl_workspace setup_ws();
gsl_workspace setup_ws(double beta);
