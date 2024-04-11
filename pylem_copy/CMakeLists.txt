cmake_minimum_required(VERSION 2.6)
project(pypf)

set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3")

find_package(PythonInterp)
find_package(PythonLibs)

find_path(BOOST_DIR "boost")

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake )

include(UseCython)
include(FindNumPy)

include_directories(../Barnes2013-Depressions/richdem/include)
include_directories("${BOOST_DIR}")
include_directories("${PYTHON_NUMPY_INCLUDE_DIR}")
include_directories("./")

set_source_files_properties(pypf.pyx PROPERTIES CYTHON_IS_CXX 1)
cython_add_module(pypf pypf.pyx pypfc.cpp)
target_link_libraries(pypf)

project(pyas)

set (CMAKE_CXX_STANDARD 11)

find_package(PythonInterp)
find_package(PythonLibs)

find_path(BOOST_DIR "boost")

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake )

include(UseCython)
include(FindNumPy)

include_directories(../Barnes2013-Depressions/richdem/include)
include_directories("${BOOST_DIR}")
include_directories("${PYTHON_NUMPY_INCLUDE_DIR}")
include_directories("./")

set_source_files_properties(pyas.pyx PROPERTIES CYTHON_IS_CXX 1)
cython_add_module(pyas pyas.pyx pyasc.cpp)
target_link_libraries(pyas)
