set(SOURCES
preparation.f90
input_Ac.f90
set_deps.f90
set_ddeps_dkz.f90
)

add_library(${PREP_LIB} STATIC ${SOURCES})
add_dependencies(${PREP_LIB} ${MOD_LIB})
target_link_libraries(${PREP_LIB} ${MOD_LIB})
