# Add compile flags and defintions to the target
# according to the options selected by the user.
#
# Usage:
# ttk_set_compile_options(<library_name>)
#
function(ttk_set_compile_options library)
	target_compile_options(${library} PRIVATE -Wall)

	if (TTK_ENABLE_KAMIKAZE)
		target_compile_definitions(${library} PUBLIC TTK_ENABLE_KAMIKAZE)
	endif()

	if (TTK_ENABLE_CPU_OPTIMIZATION)
		target_compile_options(${library} PRIVATE -march=native -O3)
	endif()

	if (TTK_ENABLE_OPENMP)
		target_compile_definitions(${library} PUBLIC TTK_ENABLE_OPENMP)
		target_compile_options(${library} PUBLIC ${OpenMP_CXX_FLAGS})
		target_link_libraries(${library} PUBLIC ${OpenMP_CXX_FLAGS})
	endif()

	if (TTK_ENABLE_MPI)
		target_compile_definitions(${library} PUBLIC TTK_ENABLE_MPI)
		target_include_directories(${library} PUBLIC ${MPI_CXX_INCLUDE_PATH})
		target_link_libraries(${library} PUBLIC ${MPI_CXX_LIBRARIES})
	endif()
endfunction()

function(ttk_print_summary)
	message(STATUS "ttk options -----------------------------------------------------------------")
	message(STATUS "TTK_ENABLE_CPU_OPTIMIZATION: ${TTK_ENABLE_CPU_OPTIMIZATION}")
	message(STATUS "TTK_ENABLE_KAMIKAZE: ${TTK_ENABLE_KAMIKAZE}")
	message(STATUS "TTK_ENABLE_OPENMP: ${TTK_ENABLE_OPENMP}")
	message(STATUS "TTK_ENABLE_MPI: ${TTK_ENABLE_MPI}")
  message(STATUS "ttk deps --------------------------------------------------------------------")
  if(Boost_INCLUDE_DIR)
    message(STATUS "Boost headers: ${Boost_INCLUDE_DIR}")
  endif()
  if(PYTHON_NUMPY_INCLUDE_DIR)
    message(STATUS "Numpy headers: ${PYTHON_NUMPY_INCLUDE_DIR}")
  endif()
  if(TTK_PYTHON_MAJOR_VERSION)
    message(STATUS 
      "Python version: ${TTK_PYTHON_MAJOR_VERSION}.${TTK_PYTHON_MINOR_VERSION}")
  endif()
	message(STATUS "ttk build -------------------------------------------------------------------")
	message(STATUS "TTK_BUILD_DOCUMENTATION: ${TTK_BUILD_DOCUMENTATION}")
	if(TTK_BUILD_DOCUMENTATION)
		message(STATUS "  DOXYGEN_EXECUTABLE: ${DOXYGEN_EXECUTABLE}")
	endif()
	message(STATUS "TTK_BUILD_STANDALONE_APPS: ${TTK_BUILD_STANDALONE_APPS}")
	message(STATUS "TTK_BUILD_VTK_WRAPPERS: ${TTK_BUILD_VTK_WRAPPERS}")
	if(TTK_BUILD_VTK_WRAPPERS)
		message(STATUS "  VTK_DIR: ${VTK_DIR}")
	endif()
	message(STATUS "ttk install -----------------------------------------------------------------")
	message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
	message(STATUS "-----------------------------------------------------------------------------")

endfunction()

function(ttk_find_python)
  find_package(PythonLibs REQUIRED QUIET)
  include_directories(${PYTHON_INCLUDE_DIRS})

  if (${CMAKE_VERSION} VERSION_GREATER "3.10" 
      OR ${CMAKE_VERSION} VERSION_EQUAL "3.10")
    string(REPLACE \".\" \" \" PYTHON_VERSION_LIST ${PYTHONLIBS_VERSION_STRING})
  else()
    string(REPLACE "." " " PYTHON_VERSION_LIST ${PYTHONLIBS_VERSION_STRING})
  endif()
  separate_arguments(PYTHON_VERSION_LIST)
  list(GET PYTHON_VERSION_LIST 0 PYTHON_MAJOR_VERSION)
  list(GET PYTHON_VERSION_LIST 1 PYTHON_MINOR_VERSION)

  set(TTK_PYTHON_MAJOR_VERSION "${PYTHON_MAJOR_VERSION}" 
    CACHE INTERNAL "TTK_PYTHON_MAJOR_VERSION")
  set(TTK_PYTHON_MINOR_VERSION "${PYTHON_MINOR_VERSION}"
    CACHE INTERNAL "TTK_PYTHON_MINOR_VERSION")

  find_path(PYTHON_NUMPY_INCLUDE_DIR numpy/arrayobject.h PATHS
    ${PYTHON_INCLUDE_DIRS}
    /usr/lib/python${PYTHON_MAJOR_VERSION}.${PYTHON_MINOR_VERSION}/site-packages/numpy/core/include/
    /usr/local/lib/python${PYTHON_MAJOR_VERSION}.${PYTHON_MINOR_VERSION}/site-packages/numpy/core/include)
  include_directories(${PYTHON_NUMPY_INCLUDE_DIR})
endfunction()

option(TTK_BUILD_DOCUMENTATION "Build doxygen developer documentation" OFF)
if(TTK_BUILD_DOCUMENTATION)
	find_package(Doxygen)
	if (DOXYGEN_FOUND)
		configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ttk.doxygen
			${CMAKE_CURRENT_BINARY_DIR}/ttk.doxygen)
		add_custom_target(doc ALL
			${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/ttk.doxygen
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			COMMENT "Generating API documentation with Doxygen" VERBATIM)
		install(DIRECTORY 
			${CMAKE_CURRENT_BINARY_DIR}/../../doc/html DESTINATION 
			${CMAKE_INSTALL_PREFIX}/share/doc/ttk)
		install(DIRECTORY
			${CMAKE_CURRENT_BINARY_DIR}/../../doc/img DESTINATION
			${CMAKE_INSTALL_PREFIX}/share/doc/ttk)
	endif()
endif()

if (NOT APPLE)
	if (MSVC)
  		option(TTK_ENABLE_OPENMP "Enable OpenMP support" FALSE)
	else()
  		option(TTK_ENABLE_OPENMP "Enable OpenMP support" TRUE)
	endif()
endif()
option(TTK_ENABLE_MPI "Enable MPI support" FALSE)

if (TTK_ENABLE_OPENMP)
	find_package(OpenMP REQUIRED)
endif()

if (TTK_ENABLE_MPI)
	find_package(MPI REQUIRED)
endif()

add_subdirectory(base)
if (TTK_BUILD_VTK_WRAPPERS)
	add_subdirectory(vtk)
endif()
