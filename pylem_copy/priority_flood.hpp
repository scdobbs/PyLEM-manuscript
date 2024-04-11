/**
  @file
  @brief Defines all the Priority-Flood algorithms described by Barnes (2014) "Priority-Flood: An Optimal Depression-Filling and Watershed-Labeling Algorithm for Digital Elevation Models".

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_priority_flood_hpp_
#define _richdem_priority_flood_hpp_
#include "Array2D.hpp"
#include "richdem/common/grid_cell.hpp"
#include <queue>
#include <limits>
#include <iostream>
#include <cstdlib> //Used for exit
using namespace richdem;

/**
  @brief  Modifies floating-point cell elevations to guarantee drainage.
  @author Richard Barnes (rbarnes@umn.edu)

    This version of Priority-Flood starts on the edges of the DEM and then
    works its way inwards using a priority queue to determine the lowest cell
    which has a path to the edge. The neighbours of this cell are added to the
    priority queue if they are higher. If they are lower, then their elevation
    is increased by a small amount to ensure that they have a drainage path and
    they are added to a "pit" queue which is used to flood pits. Cells which
    are higher than a pit being filled are added to the priority queue. In this
    way, pits are filled without incurring the expense of the priority queue.

  @param[in,out]  &elevations   A grid of cell elevations

  @pre
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM. Note that the _NoData_ value is assumed to
       be a negative number less than any actual data value.

  @post
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM.
    2. **elevations** has no landscape depressions, digital dams, or flats.
*/
template <class elev_t>
void priority_flood_epsilon(Array2D<elev_t> &elevations){
  GridCellZ_pq<elev_t> open;
  std::queue<GridCellZ<elev_t> > pit;
  //ProgressBar progress;
  //uint64_t processed_cells = 0;
  //uint64_t pitc            = 0;
  auto PitTop              = elevations.noData();
  int false_pit_cells      = 0;

  /*
  std::cerr<<"\nA Priority-Flood+Epsilon"<<std::endl;
  std::cerr<<"\nC Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024"<<std::endl;
  std::cerr<<"p Setting up boolean flood array matrix..."<<std::endl;
  */

  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);

  /*
  std::cerr<<"p Adding cells to the priority queue..."<<std::endl;
  */

  for(int x=0;x<elevations.width();x++){
    open.emplace(x,0,elevations(x,0) );
    open.emplace(x,elevations.height()-1,elevations(x,elevations.height()-1) );
    closed(x,0)=true;
    closed(x,elevations.height()-1)=true;
  }

  /*
  for(int y=1;y<elevations.height()-1;y++){
    open.emplace(0,y,elevations(0,y)  );
    open.emplace(elevations.width()-1,y,elevations(elevations.width()-1,y) );
    closed(0,y)=true;
    closed(elevations.width()-1,y)=true;
  }
  */

  /*
  std::cerr<<"p Performing Priority-Flood+Epsilon..."<<std::endl;
  progress.start( elevations.size() );
  */
  while(open.size()>0 || pit.size()>0){
    GridCellZ<elev_t> c;
    if(pit.size()>0 && open.size()>0 && open.top().z==pit.front().z){
      c=open.top();
      open.pop();
      PitTop=elevations.noData();
    } else if(pit.size()>0){
      c=pit.front();
      pit.pop();
      if(PitTop==elevations.noData())
        PitTop=elevations(c.x,c.y);
    } else {
      c=open.top();
      open.pop();
      PitTop=elevations.noData();
    }
    //processed_cells++;

    for(int n=1;n<=8;n++){
      int nx=c.x+dx[n];
      // Periodic BCs:
      nx = (nx == elevations.width()) ? 0 : (nx == -1) ? elevations.width()-1 : nx;
      int ny=c.y+dy[n];

      if(!elevations.inGrid(nx,ny)) continue;

      if(closed(nx,ny))
        continue;
      closed(nx,ny)=true;

      if(elevations(nx,ny)==elevations.noData())
        pit.push(GridCellZ<elev_t>(nx,ny,elevations.noData()));

      else if(elevations(nx,ny)<=nextafterf(c.z,std::numeric_limits<float>::infinity())){
        if(PitTop!=elevations.noData() && PitTop<elevations(nx,ny) && nextafterf(c.z,std::numeric_limits<float>::infinity())>=elevations(nx,ny))
          ++false_pit_cells;
        //++pitc;
        elevations(nx,ny)=nextafterf(c.z,std::numeric_limits<float>::infinity());
        pit.push(GridCellZ<elev_t>(nx,ny,elevations(nx,ny)));
      } else
        open.emplace(nx,ny,elevations(nx,ny));
    }
  }

}


///Priority-Flood+Epsilon is not available for integer data types
template<>
void priority_flood_epsilon(Array2D<uint8_t> &elevations){
  std::cerr<<"E Priority-Flood+Epsilon is only available for floating-point data types!"<<std::endl;
  exit(-1);
}

///Priority-Flood+Epsilon is not available for integer data types
template<>
void priority_flood_epsilon(Array2D<uint16_t> &elevations){
  std::cerr<<"E Priority-Flood+Epsilon is only available for floating-point data types!"<<std::endl;
  exit(-1);
}

///Priority-Flood+Epsilon is not available for integer data types
template<>
void priority_flood_epsilon(Array2D<int16_t> &elevations){
  std::cerr<<"E Priority-Flood+Epsilon is only available for floating-point data types!"<<std::endl;
  exit(-1);
}

///Priority-Flood+Epsilon is not available for integer data types
template<>
void priority_flood_epsilon(Array2D<uint32_t> &elevations){
  std::cerr<<"E Priority-Flood+Epsilon is only available for floating-point data types!"<<std::endl;
  exit(-1);
}

///Priority-Flood+Epsilon is not available for integer data types
template<>
void priority_flood_epsilon(Array2D<int32_t> &elevations){
  std::cerr<<"E Priority-Flood+Epsilon is only available for floating-point data types!"<<std::endl;
  exit(-1);
}


#endif
