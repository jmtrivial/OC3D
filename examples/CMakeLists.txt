

# create the binaries
ADD_EXECUTABLE(OC3D_debug OC3D_debug.cpp)
TARGET_LINK_LIBRARIES(OC3D_debug libtetmesh.a)

# recursive call
ADD_SUBDIRECTORY(Basic)
