add_subdirectory("./modules")
add_subdirectory("./main")
add_subdirectory("./preparation")
add_subdirectory("./propagation")
add_subdirectory("./analysis")
add_subdirectory("./multi_scale")

list(APPEND LIBRARIES 
  ${MOD_LIB}
  ${MAIN_LIB}
  ${PREP_LIB}
  ${PROP_LIB}
  ${ANALY_LIB}
  ${MS_LIB}
)

set(LIBRARIES ${LIBRARIES} PARENT_SCOPE)

# Executable binary
add_executable(${TARGET_NAME} main.f90)


target_link_libraries(${TARGET_NAME} ${LIBRARIES})
add_dependencies(${TARGET_NAME} ${LIBRARIES})
