/**
  @file
  @brief Defines a 2D array object with many convenient methods for working with raster data, along with several functions for checking file data types.

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_array_2d_hpp_
#define _richdem_array_2d_hpp_

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <typeinfo>
#include <stdexcept>
#include <limits>
#include <ctime>         //Used for timestamping output files
#include <unordered_set> //For printStamp
#include "richdem/common/version.hpp"
#include "richdem/common/constants.hpp"
using namespace richdem;

//These enable compression in the loadNative() and saveNative() methods
#ifdef WITH_COMPRESSION
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#endif

/**
  @brief  Class to hold and manipulate GDAL and native rasters
  @author Richard Barnes (rbarnes@umn.edu)

  Array2D manages a two-dimensional raster dataset. Passed a request to load
  such data, it peeks at the file header and can either load data on
  construction or wait until a later point. It can also offload data to disk.

  Array2D permits simple copy construction as well as templated copies, which
  transfer projections and geotransforms, but not the actual data. This is
  useful for say, create a flow directions raster which is homologous to a DEM.

  Array2D implements two addressing schemes: "xy" and "i". All methods are
  available in each scheme; users may use whichever is convenient. The xy-scheme
  accesses raster cells by their xy-coordinates. The i-scheme accesses cells by
  their address in a flat array. Internally, xy-addresses are converted to
  i-addresses. i-addressing is frequently faster because it reduces the space
  needed to store coordinates and requires no addressing mathematics; however,
  xy-addressing may be more intuitive. It is suggested to develop algorithms
  using xy-addressing and then convert them to i-addressing if additional speed
  is desired. The results of the two versions can then be compared against each
  other to verify that using i-addressing has not introduced any errors.
*/
template<class T>
class Array2D {
 public:
  std::string filename;             ///< TODO
  std::string basename;             ///< Filename without path or extension
  std::vector<double> geotransform; ///< Geotransform of the raster
  std::string projection;           ///< Projection of the raster
  std::string processing_history;   ///< List of commands previously run on this dataset

  //Using uint32_t for i-addressing allows for rasters of ~65535^2. These
  //dimensions fit easily within an int32_t xy-address.
  typedef int32_t  xy_t;            ///< xy-addressing data type
  typedef uint32_t i_t;             ///< i-addressing data type

  static const i_t NO_I = std::numeric_limits<i_t>::max();

 private:
  template<typename> friend class Array2D;

  std::vector<T> data;              ///< Holds the raster data in a 1D array
                                    ///< this improves caching versus a 2D array

  T   no_data;                      ///< NoData value of the raster
  i_t num_data_cells = NO_I;        ///< Number of cells which are not NoData

  xy_t view_width;               ///< Height of raster in cells
  xy_t view_height;              ///< Width of raster in cells

  ///@{ A rectangular subregion of a larger raster can be extracted. These
  ///   variables store the offsets of this subregion in case the subregion
  ///   needs to be saved into a raster with other subregions
  xy_t view_xoff;
  xy_t view_yoff;
  ///@}

  ///If TRUE, loadData() loads data from the cache assuming  the Native format.
  ///Otherwise, it assumes it is loading from a GDAL file.
  bool from_cache;

 public:
  Array2D(){
    view_width   = 0;
    view_height  = 0;
    view_xoff    = 0;
    view_yoff    = 0;
  }

  /**
    @brief Creates a raster of the specified dimensions

    @param[in] width   Width of the raster
    @param[in] height  Height of the raster
    @param[in] val     Initial value of all the raster's cells. Defaults to the
                       Array2D template type's default value
  */
  Array2D(xy_t width, xy_t height, const T& val = T()) : Array2D() {
    resize(width,height,val);
  }

  /**
    @brief Create a raster with the same properties and dimensions as another
           raster. No data is copied between the two.

    @param[in] other   Raster whose properties and dimensions should be copied
    @param[in] val     Initial value of all the raster's cells. Defaults to the
                       Array2D template type's default value
  */
  template<class U>
  Array2D(const Array2D<U> &other, const T& val=T()) : Array2D() {
    view_width         = other.view_width;
    view_height        = other.view_height;
    view_xoff          = other.view_xoff;
    view_yoff          = other.view_yoff;
    geotransform       = other.geotransform;
    processing_history = other.processing_history;
    projection         = other.projection;
    basename           = other.basename;
    resize(other.width(), other.height(), val);
  }

  void setCacheFilename(const std::string &filename){
    this->filename = filename;
  }

  /**
    @brief Caches the raster data and all its properties to disk. Data is then
           purged from RAM.

    @param[in] filename File to save the data to

    @post  Calls to loadData() after this will result in data being loaded from
           the cache.
  */

  ///Returns a reference to the internal data array
  T* getData() { return data.data(); }

  std::vector<T> getDataVector() {return data; }

  ///@brief Number of cells in the DEM
  i_t size() const { return view_width*view_height; }

  ///Width of the raster
  xy_t width() const { return view_width; }

  ///Height of the raster
  xy_t height() const { return view_height; }

  ///X-Offset of this subregion of whatever raster we loaded from
  xy_t viewXoff() const { return view_xoff; }

  ///Y-Offset of this subregion of whatever raster we loaded from
  xy_t viewYoff() const { return view_yoff; }

  ///Returns TRUE if no data is present in RAM
  bool empty() const { return data.empty(); }

  ///Returns the NoData value of the raster. Cells equal to this value sould
  ///generally not be used in calculations. But note that the isNoData() method
  ///is a much better choice for testing whether a cell is NoData or not.
  T noData() const { return no_data; }

  ///Finds the minimum value of the raster, ignoring NoData cells
  T min() const {
    T minval = std::numeric_limits<T>::max();
    for(auto const x: data){
      if(x==no_data)
        continue;
      minval = std::min(minval,x);
    }
    return minval;
  }

  ///Finds the maximum value of the raster, ignoring NoData cells
  T max() const {
    T maxval = std::numeric_limits<T>::min();
    for(auto const x: data){
      if(x==no_data)
        continue;
      maxval = std::max(maxval,x);
    }
    return maxval;
  }

  /**
    @brief Replace one cell value with another throughout the raster. Can
           operate on NoData cells.

    @param[in] oldval   Value to be replaced
    @param[in] newval   Value to replace 'oldval' with
  */
  void replace(const T oldval, const T newval){
    for(auto &x: data)
      if(x==oldval)
        x = newval;
  }

  /**
    @brief Counts the number of occurrences of a particular value in the raster.
           Can operate on NoData cells.

    @param[in] val   Value to be be counted

    @return The number of times 'val' appears in the raster. Will be 0 if raster
            is not loaded in RAM.
  */
  i_t countval(const T val) const {
    //TODO: Warn if raster is empty?
    i_t count=0;
    for(const auto x: data)
      if(x==val)
        count++;
    return count;
  }

  /**
    @brief Convert from index coordinates to x,y coordinates

    @param[in]  i   Index coordinate
    @param[out] x   X-coordinate of i
    @param[out] y   Y-coordinate of i
  */
  void iToxy(const i_t i, xy_t &x, xy_t &y) const {
    x = i%view_width;
    y = i/view_width;
  }

  /**
    @brief Convert from x,y coordinates to index coordinates

    @param[in]  x   X-coordinate to convert
    @param[in]  y   Y-coordinate to convert

    @return Returns the index coordinate i of (x,y)
  */
  i_t xyToI(xy_t x, xy_t y) const {
    return (i_t)y*(i_t)view_width+(i_t)x;
  }

  /**
    @brief Given a cell identified by an i-coordinate, return the i-coordinate
           of the neighbour identified by dx,dy

    @param[in]  i   i-coordinate of cell whose neighbour needs to be identified
    @param[in] dx   x-displacement of the neighbour from i
    @param[in] dy   y-displacement of the neighbour from i

    @return i-coordinate of the neighbour. Usually referred to as 'ni'
  */
  i_t nToI(i_t i, xy_t dx, xy_t dy) const {
    int32_t x=i%view_width+dx;
    int32_t y=i/view_width+dy;
    if(x<0 || y<0 || x>=view_width || y>=view_height)
      return NO_I;
    return xyToI(x,y);
  }

  /**
    @brief Given a cell identified by an i-coordinate, return the i-coordinate
           of the neighbour identified by n

    @param[in]  i   i-coordinate of cell whose neighbour needs to be identified
    @param[in]  n   Neighbour to be identified

    @return i-coordinate of the neighbour. Usually referred to as 'ni'
  */
  i_t getN(i_t i, uint8_t n) const {
    assert(0<=n && n<=8);
    xy_t x = i%view_width+(xy_t)dx[n];
    xy_t y = i/view_width+(xy_t)dy[n];
    if(x<0 || y<0 || x>=view_width || y>=view_height)
      return NO_I;
    return xyToI(x,y);
  }

  /**
    @brief Copies all the properties AND data of another raster into this one

    @param[in]  o   Raster to copy

    @return Returns this raster, with the other raster's data and properties
            copied in.
  */
  template<class U>
  T& operator=(const Array2D<U> &o){
    data               = std::vector<T>(o.data.begin(),o.data.end());
    view_height        = o.view_height;
    view_width         = o.view_width;
    view_xoff          = o.view_xoff;
    view_yoff          = o.view_yoff;
    num_data_cells     = o.num_data_cells;
    geotransform       = o.geotransform;
    projection         = o.projection;
    processing_history = o.processing_history;
    no_data            = (T)o.no_data;
    return *this;
  }

  /**
    @brief Determine if two rasters are equivalent based on dimensions,
           NoData value, and their data
  */
  bool operator==(const Array2D<T> &o){
    if(width()!=o.width() || height()!=o.height())
      return false;
    if(noData()!=o.noData())
      return false;
    return data==o.data;
  }

  /**
    @brief Whether or not a cell is NoData using x,y coordinates

    @param[in]  x   X-coordinate of cell to test
    @param[in]  y   Y-coordinate of cell to test

    @return Returns TRUE if the cell is NoData
  */
  bool isNoData(xy_t x, xy_t y) const {
    assert(0<=x && x<view_width);
    assert(0<=y && y<view_height);
    return data[xyToI(x,y)]==no_data;
  }

  /**
    @brief Whether or not a cell is NoData using i coordinates

    @param[in]  i   i-coordinate of cell to test

    @return Returns TRUE if the cell is NoData
  */
  bool isNoData(i_t i) const {
    assert(0<=i && i<size());
    return data[i]==no_data;
  }

  /**
    @brief Flips the raster from top to bottom
  */
  void flipVert(){
    for(xy_t y=0;y<view_height/2;y++)
      std::swap_ranges(
        data.begin()+xyToI(0,y),
        data.begin()+xyToI(view_width,y),
        data.begin()+xyToI(0,view_height-1-y)
      );
  }

  /**
    @brief Flips the raster from side-to-side
  */
  void flipHorz(){
    for(xy_t y=0;y<view_height;y++)
      std::reverse(data.begin()+xyToI(0,y),data.begin()+xyToI(view_width,y));
  }

  /**
    @brief Flips the raster about its diagonal axis, like a matrix tranpose.
  */
  void transpose(){
    std::cerr<<"transpose() is an experimental feature."<<std::endl;
    std::vector<T> new_data(view_width*view_height);
    for(xy_t y=0;y<view_height;y++)
    for(xy_t x=0;x<view_width;x++)
      new_data[(i_t)x*(i_t)view_height+(i_t)y] = data[xyToI(x,y)];
    data = new_data;
    std::swap(view_width,view_height);
    //TODO: Offsets?
  }

  /**
    @brief Test whether a cell lies within the boundaries of the raster

    @param[in]  x   X-coordinate of cell to test
    @param[in]  y   Y-coordinate of cell to test

    @return TRUE if cell lies within the raster
  */
  bool inGrid(xy_t x, xy_t y) const {
    return 0<=x && x<view_width && 0<=y && y<view_height;
  }

  /**
    @brief Test whether a cell lies within the boundaries of the raster.

    Obviously this bears some difference from `inGrid(x,y)`.

    @param[in]  i   i-coordinate of cell to test

    @return TRUE if cell lies within the raster
  */
  bool inGrid(i_t i) const {
    return 0<=i && i<size();
  }

  /**
    @brief Test whether a cell lies on the boundary of the raster

    @param[in]  x   X-coordinate of cell to test
    @param[in]  y   X-coordinate of cell to test

    @return TRUE if cell lies on the raster's boundary
  */
  bool isEdgeCell(xy_t x, xy_t y) const {
    return x==0 || y==0 || x==view_width-1 || y==view_height-1;
  }

  bool isTopLeft    (xy_t x, xy_t y) const { return x==0         && y==0;          }
  bool isTopRight   (xy_t x, xy_t y) const { return x==width()-1 && y==0;          }
  bool isBottomLeft (xy_t x, xy_t y) const { return x==0         && y==height()-1; }
  bool isBottomRight(xy_t x, xy_t y) const { return x==width()-1 && y==height()-1; }

  bool isTopRow    (xy_t x, xy_t y) const { return y==0;          }
  bool isBottomRow (xy_t x, xy_t y) const { return y==height()-1; }
  bool isLeftCol   (xy_t x, xy_t y) const { return x==0;          }
  bool isRightCol  (xy_t x, xy_t y) const { return x==width()-1;  }

  /**
    @brief Test whether a cell lies on the boundary of the raster

    @param[in]  i   i-coordinate of cell to test

    @return TRUE if cell lies on the raster's boundary
  */
  bool isEdgeCell(i_t i) const {
    xy_t x,y;
    iToxy(i,x,y);
    return isEdgeCell(x,y);
  }

  /**
    @brief Sets the NoData value of the raster

    @param[in]   ndval    Value to change NoData to
  */
  void setNoData(const T &ndval){
    no_data = ndval;
  }

  /**
    @brief Sets all of the raster's cells to 'val'

    @param[in]   val      Value to change the cells to
  */
  void setAll(const T val){
    std::fill(data.begin(),data.end(),val);
  }

  /**
    @brief Resize the raster. Note: this clears all the raster's data.

    @param[in]   width    New width of the raster
    @param[in]   height   New height of the raster
    @param[in]   val      Value to set all the cells to. Defaults to the
                          raster's template type default value
  */
  void resize(xy_t width, xy_t height, const T& val = T()){
    data.resize(width*height);
    setAll(val);
    view_height = height;
    view_width  = width;
  }

  /*
    @brief Resize a raster to copy another raster's dimensions. Copy properies.

    @param[in]   other    Raster to match sizes with
    @param[in]   val      Value to set all the cells to. Defaults to the
                          raster's template type default value
  */
  template<class U>
  void resize(const Array2D<U> &other, const T& val = T()){
    resize(other.width(), other.height(), val);
    geotransform       = other.geotransform;
    projection         = other.projection;
    processing_history = other.processing_history;
  }

  /**
    @brief Makes a raster larger and retains the raster's old data, similar to resize.

    Note: Using this command requires RAM equal to the sum of the old raster and
    the new raster. The old raster is placed in the upper-left of the new
    raster.

    @param[in] new_width  New width of the raster. Must be >= the old width.
    @param[in] new_height New height of the raster. Must be >= the old height.
    @param[in] val        Value to set the new cells to
  */
  void expand(xy_t new_width, xy_t new_height, const T val){
    if(new_width<view_width)
      throw std::runtime_error("expand(): new_width<view_width");
    if(new_height<view_height)
      throw std::runtime_error("expand(): new_height<view_height");

    xy_t old_width  = width();
    xy_t old_height = height();

    std::vector<T> old_data = std::move(data);

    resize(new_width,new_height,val);

    for(xy_t y=0;y<old_height;y++)
    for(xy_t x=0;x<old_width;x++)
      data[y*new_width+x] = old_data[y*old_width+x];
  }

  /**
    @brief Counts the number of cells which are not NoData.
  */
  void countDataCells(){
    num_data_cells = 0;
    for(const auto x: data)
      if(x!=no_data)
        num_data_cells++;
  }

  /**
    @brief Returns the number of cells which are not NoData. May count them.

    @return Returns the number of cells which are not NoData.
  */
  i_t numDataCells(){
    if(num_data_cells==NO_I)
      countDataCells();
    return num_data_cells;
  }

  /**
    @brief Returns the number of cells which are not NoData. Does not count them.

    countDataCells() should be call prior to running this method, or the
    non-const version of the method should be used.

    @return Returns the number of cells which are not NoData.
  */
  i_t numDataCells() const {
    return num_data_cells;
  }

  /**
    @brief Return cell value based on i-coordinate

    @param[in]   i    i-coordinate of cell whose data should be fetched.

    @return The value of the cell identified by 'i'
  */
  T& operator()(i_t i){
    assert(i>=0);
    assert(i<(i_t)view_width*view_height);
    return data[i];
  }

  /**
    @brief Return cell value based on i-coordinate

    @param[in]   i    i-coordinate of cell whose data should be fetched.

    @return The value of the cell identified by 'i'
  */
  T operator()(i_t i) const {
    assert(i>=0);
    assert(i<(i_t)view_width*view_height);
    return data[i];
  }

  /**
    @brief Return cell value based on x,y coordinates

    @param[in]   x    X-coordinate of cell whose data should be fetched.
    @param[in]   y    Y-coordinate of cell whose data should be fetched.

    @return The value of the cell identified by x,y
  */
  T& operator()(xy_t x, xy_t y){
    assert(x>=0);
    assert(y>=0);
    //std::cerr<<"Width: "<<width()<<" Height: "<<height()<<" x: "<<x<<" y: "<<y<<std::endl;
    assert(x<width());
    assert(y<height());
    return data[xyToI(x,y)];
  }

  /**
    @brief Return cell value based on x,y coordinates

    @param[in]   x    X-coordinate of cell whose data should be fetched.
    @param[in]   y    Y-coordinate of cell whose data should be fetched.

    @return The value of the cell identified by x,y
  */
  T operator()(xy_t x, xy_t y) const {
    assert(x>=0);
    assert(y>=0);
    assert(x<width());
    assert(y<height());
    return data[xyToI(x,y)];
  }

  /**
    @brief Returns a copy of the top row of the raster

    @return A vector containing a copy of the top row of the raster
  */
  std::vector<T> topRow() const {
    return getRowData(0);
  }

  /**
    @brief Returns a copy of the bottom row of the raster

    @return A vector containing a copy of the bottom row of the raster
  */
  std::vector<T> bottomRow() const {
    return getRowData(view_height-1);
  }

  /**
    @brief Returns a copy of the left column of the raster

    Top to bottom is reoriented as left to right.

    @return A vector containing a copy of the left column of the raster
  */
  std::vector<T> leftColumn() const {
    return getColData(0);
  }

  /**
    @brief Returns a copy of the right column of the raster

    Top to bottom is reoriented as left to right.

    @return A vector containing a copy of the right column of the raster
  */
  std::vector<T> rightColumn() const {
    return getColData(view_width-1);
  }

  /**
    @brief Sets an entire row of a raster to a given value.

    @param[in]   y    The row to be set
    @param[in] val    The value to set the row to
  */
  void setRow(xy_t y, const T &val){
    std::fill(data.begin()+xyToI(y,0),data.begin()+xyToI(y,view_width),val);
  }

  /**
    @brief Sets an entire column of a raster to a given value.

    @param[in]   x    The column to be set
    @param[in] val    The value to set the column to
  */
  void setCol(xy_t x, const T &val){
    for(xy_t y=0;y<view_height;y++)
      data[xyToI(x,y)] = val;
  }

  /**
    @brief Returns a copy of an arbitrary row of the raster

    @param[in]   y    The row to retrieve

    @return A vector containing a copy of the selected row
  */
  std::vector<T> getRowData(xy_t y) const {
    return std::vector<T>(data.begin()+xyToI(0,y),data.begin()+xyToI(view_width,y));
  }

  /**
    @brief Returns a copy of an arbitrary column of the raster

    @param[in]   x    The column to retrieve

    @return A vector containing a copy of the selected column
  */
  std::vector<T> getColData(xy_t x) const {
    std::vector<T> temp(view_height);
    for(xy_t y=0;y<view_height;y++)
      temp[y]=data[xyToI(x,y)];
    return temp;
  }

  ///Clears all raster data from RAM
  void clear(){
    data.clear();
    data.shrink_to_fit();
  }

  /**
    @brief Copies the geotransform, projection, and basename of another raster

    @param[in]    other    Raster to copy from
  */
  template<class U>
  void templateCopy(const Array2D<U> &other){
    geotransform       = other.geotransform;
    projection         = other.projection;
    basename           = other.basename;
    processing_history = other.processing_history;
  }

};

#endif
