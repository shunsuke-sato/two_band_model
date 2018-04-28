set(SOURCES
single_cell.f90
)

add_library(${MAIN_LIB} STATIC ${SOURCES})
add_dependencies(${MAIN_LIB} ${MOD_LIB} ${PREP_LIB} ${PROP_LIB} ${ANALY_LIB} ${MS_LIB})

