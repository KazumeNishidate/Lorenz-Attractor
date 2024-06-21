#include "rw.h"
#include "prototypes.h"

int func(double t, const double y[], double f[],
	 void *params){
  double mu = *(double *)params;
  (void)(t); /* avoid unused parameter warning */
  (void)(mu);

  f[0] = SIG*(y[1] - y[0]);
  f[1] = y[0]*(RHO - y[2]) - y[1];
  f[2] = y[0]*y[1] - BET*y[2];
  return GSL_SUCCESS;
}

int jac(double t, const double y[], double *dfdy,
     double dfdt[], void *params){
  double mu = *(double *)params; // dumy
  (void)(t); /* avoid unused parameter warning */
  (void)(mu);
  
  gsl_matrix_view dfdy_mat
    = gsl_matrix_view_array (dfdy, 3, 3);
  gsl_matrix * m = &dfdy_mat.matrix;
  gsl_matrix_set(m, 0, 0,     y[1]);
  gsl_matrix_set(m, 0, 1,     y[0]);
  gsl_matrix_set(m, 0, 2,    -y[2]);
  gsl_matrix_set(m, 1, 0, RHO-y[2]);
  gsl_matrix_set(m, 1, 1,       -1);
  gsl_matrix_set(m, 1, 2,    -y[0]);  
  gsl_matrix_set(m, 2, 0,     -SIG);
  gsl_matrix_set(m, 2, 1,      SIG);
  gsl_matrix_set(m, 2, 2,      0.0);  
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  dfdt[2] = 0.0;
  return GSL_SUCCESS;
}
