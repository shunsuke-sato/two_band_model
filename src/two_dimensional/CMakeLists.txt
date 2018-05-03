set(SOURCES
current_2d.f90
dt_evolve_2d.f90
global_variables_2d.f90
input_Ac_2d.f90
preparation_2d.f90
set_band_velocity_2d.f90
set_deps_2d.f90
single_cell_2d.f90
)

add_library(${TWO_D_LIB} STATIC ${SOURCES})
add_dependencies(${TWO_D_LIB} ${MOD_LIB})

