set(SOURCES
global_variables.f90
mpi_mod.f90
)

add_library(${MOD_LIB} STATIC ${SOURCES})


