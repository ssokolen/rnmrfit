#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void fit_1d(const double *x, const double *y, const double *knots,
            double *p, const double *lb, const double *ub, 
            int32_t n, int32_t nl, int32_t nb, int32_t np, int32_t nk,
            const double *eq, const double *iq, 
            int32_t neq, int32_t niq);

void eval_1d(const double *x, double *y, const double *knots, const double *p, 
             int32_t n, int32_t nl, int32_t nb, int32_t np, int32_t nk);

#ifdef __cplusplus
}
#endif
