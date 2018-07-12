cmake_policy(SET CMP0048 NEW) # avoid annoing warnings

project(nldct VERSION 1.0 LANGUAGES CXX)
message(STATUS "Configuring main project: " ${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 2.8.2)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/")
set(CMAKE_VERBOSE_MAKEFILE OFF)

# dependencies ################################################################

find_package(Matlab COMPONENTS MEX_COMPILER MX_LIBRARY REQUIRED)

find_package(CBLAS REQUIRED)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Link LibFFTW
find_path (FFTW_INCLUDE_DIR fftw3.h)
find_library (FFTWF_LIBRARIES NAMES fftw3f)
link_libraries (${FFTWF_LIBRARIES})
if (NOT FFTW_INCLUDE_DIR)
  message (FATAL_ERROR "FFTW3 not found.")
endif ()

add_definitions(-D_USE_MATH_DEFINES)

# add project directories #####################################################

matlab_add_mex(NAME vnldct
	nldct.h	
	matrix_funs.h
	utils/lib_image.h
	utils/lib_videot.h
	utils/utilities.h
   utils/cmd_option.h

	nldct.cpp
	matrix_funs.cpp
	utils/lib_image.cpp
	utils/lib_videot.cpp
	utils/mt19937ar.c
	utils/utilities.cpp	

	mex.cpp

	LINK_TO ${CBLAS_LIBRARIES} ${FFTWF_LIBRARIES}
)
target_include_directories(vnldct PUBLIC ${FFTW_INCLUDE_DIR})
