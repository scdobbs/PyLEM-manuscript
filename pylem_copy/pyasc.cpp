#include "pyasc.h"
#include "area_slope.hpp"
#include "priority_flood.hpp"

using namespace richdem;
using namespace std;

void pyasc_dinf(double *dem, double dx, double *a, double *s, int32_t m, int32_t n) {

  Array2D<double> elevations(n, m, 0.0);
  Array2D<double> areas(n, m, pow(dx,2));
  Array2D<double> slopes(n, m, 0.0);

  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      elevations(j,i) = dem[i*n+j];
    }
  }

  priority_flood_epsilon(elevations);
  area_slope_dinf(elevations, dx, areas, slopes);

  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      a[i*n+j] = areas(j,i);
      s[i*n+j] = slopes(j,i);
    }
  }

}

void pyasc(double *dem, double dx, double *a, double *s, int32_t m, int32_t n) {

  Array2D<double> elevations(n, m, 0.0);
  Array2D<double> areas(n, m, pow(dx,2));
  Array2D<double> slopes(n, m, 0.0);

  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      elevations(j,i) = dem[i*n+j];
    }
  }

  priority_flood_epsilon(elevations);
  area_slope(elevations, dx, areas, slopes);

  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      a[i*n+j] = areas(j,i);
      s[i*n+j] = slopes(j,i);
    }
  }

}

void pylc(double *dem, double dx, double *l, int32_t m, int32_t n) {

  Array2D<double> elevations(n, m, 0.0);
  Array2D<double> len(n, m, 0.0);

  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      elevations(j,i) = dem[i*n+j];
    }
  }

  priority_flood_epsilon(elevations);
  length_(elevations, dx, len);

  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      l[i*n+j] = len(j,i);
    }
  }

}
