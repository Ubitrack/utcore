set(UBITRACK_TEST_NAME utcore_tests)

# All cpp files for UTCORE
file(GLOB UBITRACK_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
   "${CMAKE_CURRENT_SOURCE_DIR}/*/*.cpp"
)

# All header files for UTCORE
file(GLOB UBITRACK_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "${CMAKE_CURRENT_SOURCE_DIR}/*.h"  
   "${CMAKE_CURRENT_SOURCE_DIR}/*/*.h" 
)

# UTCORE library definition
add_executable(${UBITRACK_TEST_NAME} ${UBITRACK_SRCS} ${UBITRACK_HEADERS})

target_link_libraries(${UBITRACK_TEST_NAME} ${CONAN_LIBS} utcore)

# Apply Target Properties
if(MSVC)
  if(CMAKE_CROSSCOMPILING)
    set_target_properties(${UBITRACK_TEST_NAME} PROPERTIES LINK_FLAGS "/NODEFAULTLIB:secchk")
  endif()
  set_target_properties(${UBITRACK_TEST_NAME} PROPERTIES LINK_FLAGS "/NODEFAULTLIB:libc /DEBUG")
endif()

foreach(_flag ${UBITRACK_COMPILE_FLAGS})
  set_target_properties(${UBITRACK_TEST_NAME} PROPERTIES COMPILE_FLAGS "${_flag}")
endforeach()
foreach(_flag ${UBITRACK_LINK_FLAGS})
  set_target_properties(${UBITRACK_TEST_NAME} PROPERTIES LINK_FLAGS "${_flag}")
endforeach()
foreach(_flag ${UBITRACK_LINK_FLAGS_DEBUG})
  set_target_properties(${UBITRACK_TEST_NAME} PROPERTIES LINK_FLAGS_DEBUG "${_flag}")
endforeach()
foreach(_symb ${UBITRACK_DEFINES})
  set_target_properties(${UBITRACK_TEST_NAME} PROPERTIES DEFINE_SYMBOL ${_symb})
endforeach()

# set compiler Definitions
set_target_properties(${UBITRACK_TEST_NAME} PROPERTIES COMPILE_DEFINITIONS "${UBITRACK_COMPILE_DEFINITIONS}")

set_target_properties(${UBITRACK_TEST_NAME} PROPERTIES
  DEBUG_POSTFIX "${UBITRACK_DEBUG_POSTFIX}"
  RUNTIME_OUTPUT_DIRECTORY bin
  )

# Installation of libraries, headers and documents
install(TARGETS ${UBITRACK_TEST_NAME}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )