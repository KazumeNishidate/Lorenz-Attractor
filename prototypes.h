void egg_disp(const double y[]);

int func(double t, const double y[], double f[], void *params);
int jac(double t, const double y[], double *dfdy,
	double dfdt[], void *params);
