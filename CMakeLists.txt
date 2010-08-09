PROJECT(OC3D)

cmake_minimum_required(VERSION 2.8)

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wunused -Wno-long-long -Wno-deprecated -Wnon-virtual-dtor -Wcast-align -Wchar-subscripts -Wall -Wextra -Wpointer-arith -Wformat-security -Wunused-parameter -Woverloaded-virtual -Wshadow")

# check for ITK
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

# check for popt
FIND_PATH(POPT_INCLUDE_DIR popt.h
  /usr/local/include
  /usr/include
)

FIND_LIBRARY(POPT_LIBRARIES popt
  /usr/lib
  /usr/local/lib
)

IF (POPT_LIBRARIES)
  IF (POPT_INCLUDE_DIR)
    SET(POPT_FOUND "YES")
  ENDIF (POPT_INCLUDE_DIR)
ENDIF (POPT_LIBRARIES)

# add the include directories
INCLUDE_DIRECTORIES(. src/ contrib/TetMeshLib/include/)

LINK_DIRECTORIES(${OC3D_SOURCE_DIR}/contrib/TetMeshLib/lib/)

# recursive call
ADD_SUBDIRECTORY(examples)
