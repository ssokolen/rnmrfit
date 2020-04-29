#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void fit_1d(const double *x, const double *y, 
            double *p, const double *lb, const double *ub, 
            int32_t n, int32_t nl, int32_t nb, int32_t np,
            const double *basis);

void eval_1d(const double *x, double *y, const double *p, 
             int32_t n, int32_t nl, int32_t nb, int32_t np,
             const double *basis);

#ifdef __cplusplus
}
#endif
