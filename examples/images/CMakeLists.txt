

# check for ITK
IF(ITK_FOUND AND POPT_FOUND)
  ADD_EXECUTABLE(imageCorrection imageCorrection.cpp)
  TARGET_LINK_LIBRARIES(imageCorrection ITKBasicFilters ITKCommon ITKIO ${POPT_LIBRARIES})
ENDIF(ITK_FOUND AND POPT_FOUND)
