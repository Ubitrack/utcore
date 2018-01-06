set(UBITRACK_LIBRARY_NAME utcore)


# Instrumentation support for utcore library
OPTION(ENABLE_TRACING "Enable instrumentation support?" OFF)

# Instrumentation of utcore
SET(INSTRUMENTATION_SOURCE_FOLDER "${CMAKE_SOURCE_DIR}/misc/tracing")

SET(ENABLE_TRACING_DTRACE OFF)
SET(ENABLE_TRACING_ETW OFF)
SET(ENABLE_TRACING_LTTNGUST OFF)

# Tracing Providers for various platforms
IF(APPLE AND ENABLE_TRACING)
    MESSAGE( STATUS "Enabled DTrace.")
    INCLUDE(dtrace)
    add_definitions("-DHAVE_DTRACE")
    add_definitions("-DENABLE_EVENT_TRACING")
    include_directories(${CMAKE_BINARY_DIR})
    SET(ENABLE_TRACING_DTRACE ON)
ENDIF(APPLE AND ENABLE_TRACING)

IF(WIN32 AND ENABLE_TRACING)
    MESSAGE( STATUS "Enabled Event Tracing for Windows.")
    INCLUDE(etw)
    add_definitions("-DHAVE_ETW")
    add_definitions("-DENABLE_EVENT_TRACING")
    include_directories(${CMAKE_BINARY_DIR})
    SET(ENABLE_TRACING_ETW ON)
ENDIF(WIN32 AND ENABLE_TRACING)

# example provided here:
# https://github.com/lttng/lttng-ust/tree/master/doc/examples/cmake-multiple-shared-libraries
IF(CMAKE_SYSTEM_NAME MATCHES "Linux" AND ENABLE_TRACING)
    MESSAGE( STATUS "Enabled LTTNG-UST.")
    include (FindLTTngUST REQUIRED)
    add_definitions("-DHAVE_LTTNGUST")
    add_definitions("-DENABLE_EVENT_TRACING")
    SET(ENABLE_TRACING_LTTNGUST ON)
ENDIF(CMAKE_SYSTEM_NAME MATCHES "Linux" AND ENABLE_TRACING)

set(tracing_hdr_files "")
set(tracing_src_files "")

set(tracing_extra_libraries "")

IF(ENABLE_TRACING_DTRACE)
    set(tracing_hdr_files ${tracing_hdr_files} "${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_dtrace.h")
ENDIF(ENABLE_TRACING_DTRACE)

IF(ENABLE_TRACING_ETW)
    set(tracing_hdr_files ${tracing_hdr_files} "${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_etw.h")
    set(tracing_src_files ${tracing_src_files} "${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_etw.rc")
ENDIF(ENABLE_TRACING_ETW)

IF(ENABLE_TRACING_LTTNGUST)
    include_directories(${LTTNGUST_INCLUDE_DIRS})
    set(tracing_extra_libraries LTTng::UST)
ENDIF(ENABLE_TRACING_LTTNGUST)


# Source and Header Files

# All cpp files for UTCORE
FILE(GLOB UBITRACK_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
   "${CMAKE_CURRENT_SOURCE_DIR}/*/*.cpp"
   "${CMAKE_CURRENT_SOURCE_DIR}/*/*/*.cpp"
   "${CMAKE_CURRENT_SOURCE_DIR}/*/*/*/*.cpp"
   "${CMAKE_CURRENT_SOURCE_DIR}/*/*/*/*/*.cpp"
)

# All header files for UTCORE
FILE(GLOB UBITRACK_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "${CMAKE_CURRENT_SOURCE_DIR}/*.h"  
   "${CMAKE_CURRENT_SOURCE_DIR}/*/*.h" 
   "${CMAKE_CURRENT_SOURCE_DIR}/*/*/*.h"
   "${CMAKE_CURRENT_SOURCE_DIR}/*/*/*/*.h"
   "${CMAKE_CURRENT_SOURCE_DIR}/*/*/*/*/*.h"
)

# UTCORE library definition
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(${UBITRACK_LIBRARY_NAME} SHARED ${UBITRACK_SRCS} ${tracing_src_files} ${UBITRACK_HEADERS} ${tracing_hdr_files})
ELSE(BUILD_SHARED_LIBS)
  ADD_LIBRARY(${UBITRACK_LIBRARY_NAME} STATIC ${UBITRACK_SRCS} ${tracing_src_files} ${UBITRACK_HEADERS} ${tracing_hdr_files})
ENDIF(BUILD_SHARED_LIBS)

target_link_libraries(${UBITRACK_LIBRARY_NAME} ${CONAN_LIBS} ${tracing_extra_libraries})

target_include_directories(${UBITRACK_LIBRARY_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)

if(WIN32)
  IF(BUILD_SHARED_LIBS)
    string(TOUPPER "${UBITRACK_LIBRARY_NAME}" _libname_upper)
    set(UBITRACK_COMPILE_DEFINITIONS ${UBITRACK_COMPILE_DEFINITIONS} "${_libname_upper}_DLL")
  ENDIF(BUILD_SHARED_LIBS)
endif(WIN32)

# Apply Target Properties
if(MSVC)
  if(CMAKE_CROSSCOMPILING)
    set_target_properties(${UBITRACK_LIBRARY_NAME} PROPERTIES LINK_FLAGS "/NODEFAULTLIB:secchk")
  endif()
  set_target_properties(${UBITRACK_LIBRARY_NAME} PROPERTIES LINK_FLAGS "/NODEFAULTLIB:libc /DEBUG")
endif()

foreach(_flag ${UBITRACK_COMPILE_FLAGS})
  set_target_properties(${UBITRACK_LIBRARY_NAME} PROPERTIES COMPILE_FLAGS "${_flag}")
endforeach()
foreach(_flag ${UBITRACK_LINK_FLAGS})
  set_target_properties(${UBITRACK_LIBRARY_NAME} PROPERTIES LINK_FLAGS "${_flag}")
endforeach()
foreach(_flag ${UBITRACK_LINK_FLAGS_DEBUG})
  set_target_properties(${UBITRACK_LIBRARY_NAME} PROPERTIES LINK_FLAGS_DEBUG "${_flag}")
endforeach()
foreach(_symb ${UBITRACK_DEFINES})
  set_target_properties(${UBITRACK_LIBRARY_NAME} PROPERTIES DEFINE_SYMBOL ${_symb})
endforeach()

# set compiler Definitions
set_target_properties(${UBITRACK_LIBRARY_NAME} PROPERTIES COMPILE_DEFINITIONS "${UBITRACK_COMPILE_DEFINITIONS}")

# set fPIC
set_property(TARGET ${UBITRACK_LIBRARY_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)

set_target_properties(${UBITRACK_LIBRARY_NAME} PROPERTIES
  OUTPUT_NAME "${UBITRACK_LIBRARY_NAME}${UBITRACK_DLLVERSION}"
  DEBUG_POSTFIX "${UBITRACK_DEBUG_POSTFIX}"
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_INSTALL_LIBDIR}
  )

if(NOT ANDROID)
  set_target_properties(${UBITRACK_LIBRARY_NAME} PROPERTIES
    VERSION ${UBITRACK_LIBVERSION}
    SOVERSION ${UBITRACK_SOVERSION}
    )
endif(NOT ANDROID)

# Generate tracing helper files
IF(ENABLE_TRACING_DTRACE)
    DTRACE_INSTRUMENT(${UBITRACK_LIBRARY_NAME})
ENDIF(ENABLE_TRACING_DTRACE)

IF(ENABLE_TRACING_ETW)
    ETW_INSTRUMENT(${UBITRACK_LIBRARY_NAME})
ENDIF(ENABLE_TRACING_ETW)

# Installation of libraries, headers and documents
install(TARGETS ${UBITRACK_LIBRARY_NAME}
  EXPORT library_targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

foreach(hdr ${UBITRACK_HEADERS})
  string(REGEX REPLACE "${CMAKE_BINARY_DIR}/" "" hdr2 "${hdr}")
  GET_FILENAME_COMPONENT(fpath ${hdr2} PATH)
  IF(fpath)
    install(FILES ${hdr} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${fpath}" COMPONENT dev)
  ELSE(fpath)
    install(FILES ${hdr} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT dev)
  ENDIF(fpath)
endforeach()

foreach(hdr ${tracing_hdr_files})
  string(REGEX REPLACE "${CMAKE_BINARY_DIR}/" "" hdr2 "${hdr}")
  GET_FILENAME_COMPONENT(fpath ${hdr2} PATH)
  IF(fpath)
    install(FILES ${hdr} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${fpath}" COMPONENT dev)
  ELSE(fpath)
    install(FILES ${hdr} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT dev)
  ENDIF(fpath)
endforeach()

install(EXPORT library_targets
  NAMESPACE utcore::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/utcore
  )