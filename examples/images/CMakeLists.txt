

# check for ITK
IF(ITK_FOUND)
  ADD_EXECUTABLE(imageCorrection imageCorrection.cpp)
  TARGET_LINK_LIBRARIES(imageCorrection ITKBasicFilters ITKCommon ITKIO libtetmesh.a)
ENDIF(ITK_FOUND)
