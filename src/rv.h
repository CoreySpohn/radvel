void rv_calc(double rv[], const double t[], const double params[],
             const char basis[], const int ntimes, const int nplan,
             const int specific_planet);

double timetrans_to_timeperi_c(double tc, double per, double ecc, double omega);

void model_call(double mod[], const double t[], const double params[],
                const char basis[], const double time_base, const int jit_ind,
                const int gamma_ind, const int ntimes, const int nplan,
                const int specific_planet);

void residuals(double res[], const double t[], const double rv[],
               const double rv_err[], const double params[], const char basis[],
               const double time_base, const int jit_ind, const int gamma_ind,
               const int ntimes, const int nplan);

double logprob(const double t[], const double rv[], const double rv_err[],
               const double params[], const char basis[],
               const double time_base, const int jit_ind, const int gamma_ind,
               const int ntimes, const int nplan);
