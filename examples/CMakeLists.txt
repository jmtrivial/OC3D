

# create the binaries
ADD_EXECUTABLE(ex_4genus ex_4genus.cpp)
TARGET_LINK_LIBRARIES(ex_4genus libtetmesh.a)

ADD_EXECUTABLE(OC3D_debug OC3D_debug.cpp)
TARGET_LINK_LIBRARIES(OC3D_debug libtetmesh.a)