set(SOURCES
global_variables_ms.f90
preparation_ms.f90
current_ms.f90
input_Ac_ms.f90
prop_Ac_ms.f90
prop_elec_ms.f90
dt_evolve_ms.f90
multi_scale.f90
write_field_iter.f90
write_results.f90
energy_calc.f90
)

add_library(${MS_LIB} STATIC ${SOURCES})
add_dependencies(${MS_LIB} ${MOD_LIB})
target_link_libraries(${MS_LIB} ${MOD_LIB})
