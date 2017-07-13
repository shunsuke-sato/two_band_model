cmake_minimum_required(VERSION 2.8)

project(PANDA Fortran)

set(CMAKE_Fortran_COMPILER      "mpif90")
set(CMAKE_Fortran_FLAGS   "-O0 ")

set(TARGET_NAME   "panda")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

set(PANDA_SOURCES   src/main.f90 src/single_cell.f90)
file(GLOB MOD_SRC      src/modules/*.f90)
file(GLOB PREP_SRC     src/preparation/*.f90)
file(GLOB PROP_SRC     src/propagation/*.f90)
file(GLOB ANA_SRC      src/analysis/*.f90)
file(GLOB MS_SRC       src/multi_scale/*.f90)
list(APPEND PANDA_SOURCES 
            ${MOD_SRC}
            ${PREP_SRC}
            ${PROP_SRC}
            ${ANA_SRC}
            ${MS_SRC}
)

add_executable(${TARGET_NAME} ${PANDA_SOURCES})
