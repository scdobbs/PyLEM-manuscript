#include "pypfc.h"
#include "priority_flood.hpp"

using namespace std;

void pypfc(double *dem, int32_t m, int32_t n) {

  Array2D<double> elevations(n, m, 0.0);

  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      elevations(j,i) = dem[i*n+j];
    }
  }

  priority_flood_epsilon(elevations);

  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      dem[i*n+j] = elevations(j,i);
    }
  }

}
