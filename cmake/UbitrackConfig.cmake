# Version of UBITRACK API.
SET( UBITRACK_MAJOR_VERSION 1 )
SET( UBITRACK_MINOR_VERSION 3 )
SET( UBITRACK_PATCH_VERSION 0 )
SET( UBITRACK_FULL_VERSION "${UBITRACK_MAJOR_VERSION}.${UBITRACK_MINOR_VERSION}.${UBITRACK_PATCH_VERSION}")


# Target Properties
set(UBITRACK_COMPILE_FLAGS "")
set(UBITRACK_COMPILE_DEFINITIONS "")
set(UBITRACK_LINK_FLAGS "")
set(UBITRACK_LINK_FLAGS_DEBUG "")

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# do we require C++11 compilers now ?
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set(UBITRACK_COMPILE_DEFINITIONS ${UBITRACK_COMPILE_DEFINITIONS} COMPILER_USE_CXX11)
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please consider upgrading your C++ compiler.")
else()
    message(WARNING "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# need to review these settings if they are still appropriate
MESSAGE(STATUS "Building for ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}")
if(WIN32)
  set(UBITRACK_COMPILE_FLAGS ${UBITRACK_COMPILE_FLAGS} "/EHsc /c /W3 /GR /wd4355 /wd4996 /wd4251 /wd4275 /wd4819 /wd4290")
  set(UBITRACK_LINK_FLAGS ${UBITRACK_LINK_FLAGS} "/SUBSYSTEM:CONSOLE")
  set(UBITRACK_COMPILE_DEFINITIONS ${UBITRACK_COMPILE_DEFINITIONS} "WIN32" "_MBCS" "BOOST_SPIRIT_USE_OLD_NAMESPACE")  
  set(UBITRACK_LINK_FLAGS_DEBUG ${UBITRACK_LINK_FLAGS_DEBUG} "/NODEFAULTLIB:libc.lib /NODEFAULTLIB:libcmt.lib /NODEFAULTLIB:msvcrt.lib /NODEFAULTLIB:libcd.lib /NODEFAULTLIB:libcmtd.lib")

  ## Check for Windows Version ##
  if( ${CMAKE_SYSTEM_VERSION} EQUAL 6.1 ) # Windows 7
    MESSAGE(STATUS "Setting minimum Windows version to Win7 WINVER=0x0601")
     set(UBITRACK_COMPILE_DEFINITIONS ${UBITRACK_COMPILE_DEFINITIONS} WINVER=0x0601)
  elseif( ${CMAKE_SYSTEM_VERSION} EQUAL 6.2 ) # Windows 8
    MESSAGE(STATUS "Setting minimum Windows version to Win8 WINVER=0x0602")
    set(UBITRACK_COMPILE_DEFINITIONS ${UBITRACK_COMPILE_DEFINITIONS} WINVER=0x0602)
  elseif( ${CMAKE_SYSTEM_VERSION} EQUAL 6.3 ) # Windows 8.1
    MESSAGE(STATUS "Setting minimum Windows version to Win8.1 WINVER=0x0603")
    set(UBITRACK_COMPILE_DEFINITIONS ${UBITRACK_COMPILE_DEFINITIONS} WINVER=0x0603)
  elseif( ${CMAKE_SYSTEM_VERSION} EQUAL 10.0 ) # Windows 10
    MESSAGE(STATUS "Setting minimum Windows version to Win8.1 WINVER=0x0603")
    set(UBITRACK_COMPILE_DEFINITIONS ${UBITRACK_COMPILE_DEFINITIONS} WINVER=0x0603)
  else() # Some other Windows
    MESSAGE(STATUS "Setting minimum Windows version to Vista WINVER=0x0600")
    set(UBITRACK_COMPILE_DEFINITIONS ${UBITRACK_COMPILE_DEFINITIONS} WINVER=0x0600)
  endif()
endif(WIN32)

IF(NOT WIN32)
  MESSAGE(STATUS "Set boost::ublas alignment to 16")
  set(UBITRACK_COMPILE_DEFINITIONS ${UBITRACK_COMPILE_DEFINITIONS} "BOOST_UBLAS_BOUNDED_ARRAY_ALIGN=__attribute__ ((aligned (16)))")
ENDIF(NOT WIN32)


if(WIN32)
  # Postfix of DLLs:
  set(UBITRACK_DLLVERSION "${UBITRACK_MAJOR_VERSION}${UBITRACK_MINOR_VERSION}${UBITRACK_PATCH_VERSION}")
  set(UBITRACK_DEBUG_POSTFIX d)
else()
  # Postfix of so's:
  set(UBITRACK_DLLVERSION "")
  set(UBITRACK_DEBUG_POSTFIX "")
endif()

if(DEFINED CMAKE_DEBUG_POSTFIX)
  set(UBITRACK_DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")
endif()

# Version of project.
SET(UBITRACK_SOVERSION "${UBITRACK_MAJOR_VERSION}.${UBITRACK_MINOR_VERSION}")
SET(UBITRACK_LIBVERSION "${UBITRACK_MAJOR_VERSION}.${UBITRACK_MINOR_VERSION}.${UBITRACK_PATCH_VERSION}")

if (APPLE)
    set(CMAKE_INSTALL_RPATH "@executable_path/../lib")
else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
endif()

# default install directories
set(UBITRACK_DOC_INSTALL_DIRECTORY "${CMAKE_INSTALL_DATADIR}/UbiTrack/doc")
set(UBITRACK_UTQLPATTERN_INSTALL_DIRECTORY "${CMAKE_INSTALL_DATADIR}/Ubitrack/utql")

set(UBITRACK_COMPONENT_INSTALL_DIRECTORY "ubitrack")
set(UBITRACK_COMPONENT_INSTALL_PATH "${CMAKE_INSTALL_LIBDIR}/${UBITRACK_COMPONENT_INSTALL_DIRECTORY}")
set(UBITRACK_COMPONENT_BIN_INSTALL_PATH "${CMAKE_INSTALL_BINDIR}/${UBITRACK_COMPONENT_INSTALL_DIRECTORY}")


# Helper Marcros to simpilify and unify the installation process of libraries and components 

# Macro to set target_properties from the custom variables defined above
macro(ubitrack_set_target_properties target_name)

  foreach(_flag ${UBITRACK_COMPILE_FLAGS})
    set_target_properties(${target_name} PROPERTIES COMPILE_FLAGS "${_flag}")
  endforeach()
  foreach(_flag ${UBITRACK_LINK_FLAGS})
    set_target_properties(${target_name} PROPERTIES LINK_FLAGS "${_flag}")
  endforeach()
  foreach(_flag ${UBITRACK_LINK_FLAGS_DEBUG})
    set_target_properties(${target_name} PROPERTIES LINK_FLAGS_DEBUG "${_flag}")
  endforeach()
  foreach(_symb ${UBITRACK_DEFINES})
    set_target_properties(${target_name} PROPERTIES DEFINE_SYMBOL ${_symb})
  endforeach()

  # set compiler Definitions
  set_target_properties(${target_name} PROPERTIES COMPILE_DEFINITIONS "${UBITRACK_COMPILE_DEFINITIONS}")

  # set fPIC
  set_target_properties(${target_name} PROPERTIES POSITION_INDEPENDENT_CODE ON)

  set_target_properties(${target_name} PROPERTIES
    OUTPUT_NAME "${target_name}${UBITRACK_DLLVERSION}"
    DEBUG_POSTFIX "${UBITRACK_DEBUG_POSTFIX}"
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_INSTALL_LIBDIR}
    )

  if(NOT ANDROID)
    set_target_properties(${target_name} PROPERTIES
      VERSION ${UBITRACK_LIBVERSION}
      SOVERSION ${UBITRACK_SOVERSION}
      )
  endif(NOT ANDROID)

endmacro(ubitrack_set_target_properties)

# Macro to install the header files in a standarts conform way
# Multiple arguments are allowed: each one is a list of header files to copy
macro(ubitrack_install_headers)
  foreach(arg ${ARGN})
    foreach(hdr ${arg})
      string(REGEX REPLACE "${CMAKE_BINARY_DIR}/" "" hdr2 "${hdr}")
      GET_FILENAME_COMPONENT(fpath ${hdr2} PATH)
      IF(fpath)
        install(FILES ${hdr} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${fpath}" COMPONENT dev)
      ELSE(fpath)
        install(FILES ${hdr} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT dev)
      ENDIF(fpath)
    endforeach()
  endforeach()
endmacro(ubitrack_install_headers)

# Macro to install the current library in a standarts conform way
# One argument is allowed: target_name
macro(ubitrack_install_library target_name)
  install(TARGETS ${target_name}
    EXPORT "${target_name}Targets"
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

  write_basic_package_version_file("${CMAKE_SOURCE_DIR}/${target_name}ConfigVersion.cmake"
    VERSION ${UBITRACK_FULL_VERSION}
    COMPATIBILITY SameMajorVersion)

  # this requires the ${target_name}Config.cmake file to be present in all libraries
  install (FILES "${CMAKE_SOURCE_DIR}/${target_name}Config.cmake" "${CMAKE_SOURCE_DIR}/${target_name}ConfigVersion.cmake"
    DESTINATION lib/cmake/${target_name})

  install(EXPORT "${target_name}Targets"
    FILE "${target_name}Targets.cmake"
    NAMESPACE "${target_name}::"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${target_name}"
  )
endmacro(ubitrack_install_library)

# Macro to install the current library in a standarts conform way
# One argument is allowed: target_name
macro(ubitrack_install_component target_name)
  install(TARGETS ${target_name}
    RUNTIME DESTINATION ${UBITRACK_COMPONENT_BIN_INSTALL_PATH}
    LIBRARY DESTINATION ${UBITRACK_COMPONENT_INSTALL_PATH}
    ARCHIVE DESTINATION ${UBITRACK_COMPONENT_INSTALL_PATH}
    )
endmacro(ubitrack_install_component)

# Macro to install the current application in a standarts conform way
# One argument is allowed: target_name
macro(ubitrack_install_application target_name)
  install(TARGETS ${target_name}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
endmacro(ubitrack_install_application)

# Macro to install additional data files into the doc folder
# multiple arguments allowed: filenames or glob patterns
macro(ubitrack_install_doc)
  foreach(arg ${ARGN})
    file(GLOB _doc_files LIST_DIRECTORIES false ${arg})
    foreach(pfile ${_doc_files})
      string(REGEX REPLACE "^.*/doc/" "" pfile2 "${pfile}")
      GET_FILENAME_COMPONENT(fpath ${pfile2} PATH)
      IF(fpath)
        install(FILES ${pfile} DESTINATION "${UBITRACK_DOC_INSTALL_DIRECTORY}/${fpath}" COMPONENT doc)
      ELSE(fpath)
        install(FILES ${pfile} DESTINATION "${UBITRACK_DOC_INSTALL_DIRECTORY}" COMPONENT doc)
      ENDIF(fpath)
    endforeach()
  endforeach()
endmacro(ubitrack_install_doc)

# Macro to install utql data files into the doc folder
# multiple arguments allowed: filenames or glob patterns
macro(ubitrack_install_utql)
  foreach(arg ${ARGN})
    file(GLOB _utql_files LIST_DIRECTORIES false ${arg})
    foreach(pfile ${_utql_files})
      string(REGEX REPLACE "^.*/doc/utql/" "" pfile2 "${pfile}")
      GET_FILENAME_COMPONENT(fpath ${pfile2} PATH)
      IF(fpath)
        install(FILES ${pfile} DESTINATION "${UBITRACK_UTQLPATTERN_INSTALL_DIRECTORY}/${fpath}" COMPONENT doc)
      ELSE(fpath)
        install(FILES ${pfile} DESTINATION "${UBITRACK_UTQLPATTERN_INSTALL_DIRECTORY}" COMPONENT doc)
      ENDIF(fpath)
    endforeach()
  endforeach()
endmacro(ubitrack_install_utql)