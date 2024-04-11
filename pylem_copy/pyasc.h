#include <stdint.h>

#ifndef PYAS_H
#define PYAS_H

void pyasc(double *dem, double dx, double *a, double *s, int32_t m, int32_t n);
void pyasc_dinf(double *dem, double dx, double *a, double *s, int32_t m, int32_t n);
void pylc(double *dem, double dx, double *l, int32_t m, int32_t n);

#endif // PYPF_H
