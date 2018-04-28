set(SOURCES
energy.f90
current.f90
excited_electron.f90
excited_electron_k_resolved.f90
)

add_library(${ANALY_LIB} STATIC ${SOURCES})
add_dependencies(${ANALY_LIB} ${MOD_LIB} ${PREP_LIB} ${PROP_LIB})
