cmake_minimum_required(VERSION 2.8)



#set(CMAKE_Fortran_COMPILER      "mpif90")
set(CMAKE_Fortran_COMPILER      "mpiifort") 

#set(CMAKE_Fortran_FLAGS   "-O0 -fopenmp")
set(CMAKE_Fortran_FLAGS   "-O3 -ipo -xHOST -qopenmp")

project(PANDA Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})

set(TARGET_NAME   "panda")

set(MOD_LIB "module_lib")
set(MAIN_LIB "main_lib")
set(PREP_LIB   "preparation_lib")
set(PROP_LIB   "propagation_lib")
set(ANALY_LIB   "analysis_lib")
set(MS_LIB   "multi_scale_lib")
set(TWO_D_LIB   "two_dimensional_lib")

set(SRC_DIR "./src")
add_subdirectory(${SRC_DIR})

