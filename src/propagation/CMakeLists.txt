set(SOURCES
dt_evolve.f90
)

add_library(${PROP_LIB} STATIC ${SOURCES})
add_dependencies(${PROP_LIB} ${MOD_LIB} ${PREP_LIB})
target_link_libraries(${PROP_LIB} ${MOD_LIB} ${PREP_LIB})
