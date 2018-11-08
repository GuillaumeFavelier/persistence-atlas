cmake_minimum_required(VERSION 3.2)

project(persistenceAtlasCmd)

set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include_directories(${PYTHON_NUMPY_INCLUDE_DIR})
add_executable(persistenceAtlasCmd main.cpp)
target_link_libraries(persistenceAtlasCmd
  ttkSortScalarField
  ttkPersistenceMap
  ttkLDistanceMatrix
  ttkSpectralEmbedding
  ttkClustering
  ttkSphereFromPoint
  ttkUncertainDataEstimator
  ttkMandatoryCriticalPoints
  ttkProgramBase
  ${PYTHON_LIBRARY}
  ${VTK_LIBRARIES})

install(TARGETS persistenceAtlasCmd RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
