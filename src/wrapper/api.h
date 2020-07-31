#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void fit_1d(const double *x, const double *y, const double *knots,
            double *p, const double *lb, const double *ub, 
            int32_t n, int32_t nl, int32_t nb, int32_t np, int32_t nk,
            const double *eq, const double *iq, 
            int32_t neq, int32_t niq, int32_t alg, double xtr, double mxt);

void eval_1d(const double *x, double *y, const double *knots, const double *p, 
             int32_t n, int32_t nl, int32_t nb, int32_t np, int32_t nk);

void baseline_1d(const double *x, double *y, const double *knots, const double *p, 
                 int32_t n, int32_t nb, int32_t nk);

void phase_1d(const double *x, double *y, const double *p, int32_t n, int32_t np);

void fit_2d(const double *x_direct, const double *x_indirect, const double *y, 
            const int32_t *resonances, const int32_t *dimensions,
            const double *knots, double *p, const double *lb, const double *ub, 
            int32_t n, int32_t nl, int32_t nb, int32_t np, int32_t nk,
            const double *eq, const double *iq, 
            int32_t neq, int32_t niq);

void eval_2d(const double *x_direct, const double *x_indirect, double *y, 
             const int32_t *resonances, const int32_t *dimensions,
             const double *knots, const double *p, 
             int32_t n, int32_t nl, int32_t nb, int32_t np, int32_t nk);

void baseline_2d(const double *x_direct, const double *x_indirect, double *y, 
                 const double *knots, const double *p, 
                 int32_t n, int32_t nb, int32_t nk);

void phase_2d(const double *x_direct, const double *x_indirect, double *y, 
              const double *p, int32_t n, int32_t np);

#ifdef __cplusplus
}
#endif
