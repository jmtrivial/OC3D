PROJECT(OC3D)

cmake_minimum_required(VERSION 2.8)

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wunused -Wno-long-long -Wno-deprecated -Wnon-virtual-dtor -Wcast-align -Wchar-subscripts -Wall -Wextra -Wpointer-arith -Wformat-security -Wunused-parameter -Woverloaded-virtual -Wshadow")


# add the include directories
INCLUDE_DIRECTORIES(. src/ contrib/TetMeshLib/include/)

LINK_DIRECTORIES(${OC3D_SOURCE_DIR}/contrib/TetMeshLib/lib/)

# recursive call
ADD_SUBDIRECTORY(examples)
