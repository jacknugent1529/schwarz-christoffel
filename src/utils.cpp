#include "utils.hpp"

void free_ws(gsl_workspace ws) {
    gsl_integration_workspace_free(ws.ws);
    if (ws.table != NULL) {
        gsl_integration_qaws_table_free(ws.table);
    }
}

gsl_workspace setup_ws() {
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(100000);
    return { w, NULL };
}

gsl_workspace setup_ws(double beta) {
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(100000);
    gsl_integration_qaws_table *table = gsl_integration_qaws_table_alloc(0, beta, 0, 0);
    return { w, table};
}
