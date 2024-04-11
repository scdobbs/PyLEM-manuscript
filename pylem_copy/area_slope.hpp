#ifndef _area_slope_hpp_
#define _area_slope_hpp_
#include "Array2D.hpp"
#include "richdem/common/grid_cell.hpp"

#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include <cmath>

using namespace std;

typedef int32_t  xy_t;

template <class elev_t>
void update(xy_t this_x, xy_t this_y, xy_t next_x, xy_t next_y, elev_t dx, Array2D<elev_t> &elevations, double &maxSlope, xy_t &max_next_x, xy_t &max_next_y) {
  elev_t thisSlope = (elevations(this_x, this_y) - elevations(next_x, next_y)) / (((next_x != this_x) || (next_y != this_y)) ? (1.41*dx): dx);
  if(thisSlope > maxSlope) {
    maxSlope = thisSlope;
    max_next_x = next_x;
    max_next_y = next_y;
  }
}


template <class elev_t>
void area_slope(Array2D<elev_t> &elevations, elev_t dx, Array2D<elev_t> &area, Array2D<elev_t> &slope) {

  vector<elev_t> v = elevations.getDataVector();
  vector<size_t> indices(v.size());
  iota(indices.begin(), indices.end(), 0);
  sort(indices.begin(), indices.end(),
              [&v](int left, int right) -> bool {
                  // sort indices according to corresponding array element
                  return v[left] < v[right];
              });
  reverse(indices.begin(), indices.end());

  xy_t nx, ny;

  nx = elevations.width();
  ny = elevations.height();

  for (auto i: indices) {

    xy_t this_x, this_y, next_x, next_y, max_next_x, max_next_y;

    elevations.iToxy(i, this_x, this_y);

    elev_t maxSlope = 0;

    if( (this_y != 0) && (this_y != (ny-1))) {

      next_x = (this_x-1 == -1) ? nx - 1 : this_x - 1;
      next_y = this_y + 1;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      next_x = this_x;
      next_y = this_y + 1;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      next_x = (this_x+1 == nx) ? 0 : this_x + 1;
      next_y = this_y + 1;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      next_x = (this_x+1 == nx) ? 0 : this_x + 1;
      next_y = this_y;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      next_x = (this_x+1 == nx) ? 0 : this_x + 1;
      next_y = this_y - 1;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      next_x = this_x;
      next_y = this_y - 1;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      next_x = (this_x - 1 == -1) ? nx-1 : this_x - 1;
      next_y = this_y - 1;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      next_x = (this_x - 1 == -1) ? nx-1 : this_x - 1;
      next_y = this_y;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      area(max_next_x, max_next_y) += area(this_x, this_y);
      slope(this_x, this_y) = maxSlope;

    }

  }

}

template <class elev_t>
void update_dinf(xy_t this_x, xy_t this_y, xy_t next_x1, xy_t next_y1, xy_t next_x2, xy_t next_y2, elev_t dx,
  Array2D<elev_t> &elevations, double &maxSlope, xy_t &max_next_x1, xy_t &max_next_y1, xy_t &max_next_x2,
  xy_t &max_next_y2, elev_t &partition1, elev_t &partition2) {

  bool x1IsDiagonal = (next_x1 != this_x) && (next_y1 != this_y);
  elev_t s1, s2, r, thisSlope;

  if(!x1IsDiagonal) {
    s1 = (elevations(this_x, this_y) - elevations(next_x1, next_y1)) / dx;
    s2 = (elevations(next_x1, next_y1) - elevations(next_x2, next_y2)) / dx;
  } else {
    s1 = (elevations(this_x, this_y) - elevations(next_x2, next_y2)) / dx;
    s2 = (elevations(next_x2, next_y2) - elevations(next_x1, next_y1)) / dx;
  }

  r = atan2(s2, s1);
  if(r < 0) {
    r = 0;
    thisSlope = s1;
  } else if(r > atan2(1,1)) {
    r = atan2(1,1);
    if(x1IsDiagonal) {
      thisSlope =(elevations(this_x, this_y) - elevations(next_x1, next_y1)) / (sqrt(2)*dx);
    } else {
      thisSlope =(elevations(this_x, this_y) - elevations(next_x2, next_y2)) / (sqrt(2)*dx);
    }
  } else {
    thisSlope = sqrt(pow(s1,2) + pow(s2,2));
  }

  if(thisSlope > maxSlope) {
    maxSlope = thisSlope;
    if(!x1IsDiagonal) {
      partition1 = 1 - tan(r);
      partition2 = tan(r);
    } else {
      partition2 = 1 - tan(r);
      partition1 = tan(r);
    }
    max_next_x1 = next_x1;
    max_next_y1 = next_y1;
    max_next_x2 = next_x2;
    max_next_y2 = next_y2;
  }
}

template <class elev_t>
void area_slope_dinf(Array2D<elev_t> &elevations, elev_t dx, Array2D<elev_t> &area, Array2D<elev_t> &slope) {

  vector<elev_t> v = elevations.getDataVector();
  vector<size_t> indices(v.size());
  iota(indices.begin(), indices.end(), 0);
  sort(indices.begin(), indices.end(),
              [&v](int left, int right) -> bool {
                  // sort indices according to corresponding array element
                  return v[left] < v[right];
              });
  reverse(indices.begin(), indices.end());

  xy_t nx, ny;

  nx = elevations.width();
  ny = elevations.height();

  for (auto i: indices) {

    xy_t this_x, this_y, next_x1, next_y1, next_x2, next_y2, max_next_x1, max_next_y1, max_next_x2, max_next_y2;

    elevations.iToxy(i, this_x, this_y);

    elev_t maxSlope = -1.0;
    elev_t partition1 = 0;
    elev_t partition2 = 0;

    if( (this_y != 0) && (this_y != (ny-1))) {

      // Facet 6:

      next_x1 = (this_x-1 == -1) ? nx - 1 : this_x - 1;
      next_y1 = this_y + 1;
      next_x2 = this_x;
      next_y2 = this_y + 1;
      update_dinf(this_x, this_y, next_x1, next_y1, next_x2, next_y2, dx, elevations, maxSlope, max_next_x1, max_next_y1, max_next_x2, max_next_y2, partition1, partition2);

      // Facet 7:

      next_x1 = next_x2;
      next_y1 = next_y2;
      next_x2 = (this_x+1 == nx) ? 0 : this_x + 1;
      next_y2 = this_y + 1;
      update_dinf(this_x, this_y, next_x1, next_y1, next_x2, next_y2, dx, elevations, maxSlope, max_next_x1, max_next_y1, max_next_x2, max_next_y2, partition1, partition2);

      // Facet 8:

      next_x1 = next_x2;
      next_y1 = next_y2;
      next_x2 = (this_x+1 == nx) ? 0 : this_x + 1;
      next_y2 = this_y;
      update_dinf(this_x, this_y, next_x1, next_y1, next_x2, next_y2, dx, elevations, maxSlope, max_next_x1, max_next_y1, max_next_x2, max_next_y2, partition1, partition2);

      // Facet 1:

      next_x1 = next_x2;
      next_y1 = next_y2;
      next_x2 = (this_x+1 == nx) ? 0 : this_x + 1;
      next_y2 = this_y - 1;
      update_dinf(this_x, this_y, next_x1, next_y1, next_x2, next_y2, dx, elevations, maxSlope, max_next_x1, max_next_y1, max_next_x2, max_next_y2, partition1, partition2);

      // Facet 2:

      next_x1 = next_x2;
      next_y1 = next_y2;
      next_x2 = this_x;
      next_y2 = this_y - 1;
      update_dinf(this_x, this_y, next_x1, next_y1, next_x2, next_y2, dx, elevations, maxSlope, max_next_x1, max_next_y1, max_next_x2, max_next_y2, partition1, partition2);

      // Facet 3:

      next_x1 = next_x2;
      next_y1 = next_y2;
      next_x2 = (this_x-1 == -1) ? nx - 1 : this_x - 1;
      next_y2 = this_y - 1;
      update_dinf(this_x, this_y, next_x1, next_y1, next_x2, next_y2, dx, elevations, maxSlope, max_next_x1, max_next_y1, max_next_x2, max_next_y2, partition1, partition2);

      // Facet 4:

      next_x1 = next_x2;
      next_y1 = next_y2;
      next_x2 = (this_x-1 == -1) ? nx - 1 : this_x - 1;
      next_y2 = this_y;
      update_dinf(this_x, this_y, next_x1, next_y1, next_x2, next_y2, dx, elevations, maxSlope, max_next_x1, max_next_y1, max_next_x2, max_next_y2, partition1, partition2);

      // Facet 5:

      next_x1 = next_x2;
      next_y1 = next_y2;
      next_x2 = (this_x-1 == -1) ? nx - 1 : this_x - 1;;
      next_y2 = this_y + 1;
      update_dinf(this_x, this_y, next_x1, next_y1, next_x2, next_y2, dx, elevations, maxSlope, max_next_x1, max_next_y1, max_next_x2, max_next_y2, partition1, partition2);


      if(maxSlope > 0) {
        area(max_next_x1, max_next_y1) += area(this_x, this_y)*partition1;
        area(max_next_x2, max_next_y2) += area(this_x, this_y)*partition2;
        slope(this_x, this_y) = maxSlope;
      }
    }

  }

}


template <class elev_t>
void length_(Array2D<elev_t> &elevations, elev_t dx, Array2D<elev_t> &length) {

  vector<elev_t> v = elevations.getDataVector();
  vector<size_t> indices(v.size());
  iota(indices.begin(), indices.end(), 0);
  sort(indices.begin(), indices.end(),
              [&v](int left, int right) -> bool {
                  // sort indices according to corresponding array element
                  return v[left] < v[right];
              });
  reverse(indices.begin(), indices.end());

  xy_t nx, ny;

  nx = elevations.width();
  ny = elevations.height();

  for (auto i: indices) {

    xy_t this_x, this_y, next_x, next_y, max_next_x, max_next_y;

    elevations.iToxy(i, this_x, this_y);

    elev_t maxSlope = 0;

    if( (this_y != 0) && (this_y != (ny-1))) {

      next_x = (this_x-1 == -1) ? nx - 1 : this_x - 1;
      next_y = this_y + 1;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      next_x = this_x;
      next_y = this_y + 1;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      next_x = (this_x+1 == nx) ? 0 : this_x + 1;
      next_y = this_y + 1;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      next_x = (this_x+1 == nx) ? 0 : this_x + 1;
      next_y = this_y;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      next_x = (this_x+1 == nx) ? 0 : this_x + 1;
      next_y = this_y - 1;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      next_x = this_x;
      next_y = this_y - 1;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      next_x = (this_x - 1 == -1) ? nx-1 : this_x - 1;
      next_y = this_y - 1;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      next_x = (this_x - 1 == -1) ? nx-1 : this_x - 1;
      next_y = this_y;
      update(this_x, this_y, next_x, next_y, dx, elevations, maxSlope, max_next_x, max_next_y);

      if(length(max_next_x, max_next_y) < length(this_x, this_y) + (((this_x == next_x) || (this_y == next_y)) ? dx : 1.414*dx)) {
        length(max_next_x, max_next_y) = length(this_x, this_y) + (((this_x == next_x) || (this_y == next_y)) ? dx : 1.414*dx);
      }
    }

  }

}


#endif
