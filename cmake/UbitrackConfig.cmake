# Version of UBITRACK API.
SET( UBITRACK_MAJOR_VERSION 1 )
SET( UBITRACK_MINOR_VERSION 3 )
SET( UBITRACK_PATCH_VERSION 0 )

# do we require C++11 compilers now ?
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# Target Properties
set(UBITRACK_COMPILE_FLAGS "")
set(UBITRACK_COMPILE_DEFINITIONS "")
set(UBITRACK_LINK_FLAGS "")
set(UBITRACK_LINK_FLAGS_DEBUG "")

# need to review these settings if they are still appropriate
MESSAGE(STATUS "Building for ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}")
if(WIN32)
  set(UBITRACK_COMPILE_FLAGS "/EHsc /c /W3 /GR /wd4355 /wd4996 /wd4251 /wd4275 /wd4819 /wd4290")
  set(UBITRACK_LINK_FLAGS "/SUBSYSTEM:CONSOLE")
  set(UBITRACK_COMPILE_DEFINITIONS "WIN32" "_MBCS" "BOOST_SPIRIT_USE_OLD_NAMESPACE")  
  set(UBITRACK_LINK_FLAGS_DEBUG "/NODEFAULTLIB:libc.lib /NODEFAULTLIB:libcmt.lib /NODEFAULTLIB:msvcrt.lib /NODEFAULTLIB:libcd.lib /NODEFAULTLIB:libcmtd.lib")

  ## Check for Windows Version ##
  if( ${CMAKE_SYSTEM_VERSION} EQUAL 6.1 ) # Windows 7
    MESSAGE(STATUS "Setting minimum Windows version to Win7 WINVER=0x0601")
     set(UBITRACK_COMPILE_DEFINITIONS ${UBITRACK_COMPILE_DEFINITIONS} WINVER=0x0601)
  elseif( ${CMAKE_SYSTEM_VERSION} EQUAL 6.2 ) # Windows 8
    MESSAGE(STATUS "Setting minimum Windows version to Win8 WINVER=0x0602")
    set(UBITRACK_COMPILE_DEFINITIONS ${UBITRACK_COMPILE_DEFINITIONS} WINVER=0x0602)
  elseif( ${CMAKE_SYSTEM_VERSION} EQUAL 6.3 ) # Windows 8.1
    MESSAGE(STATUS "Setting minimum Windows version to Win8.1 WINVER=0x0603")
    set(UBITRACK_COMPILE_DEFINITIONS ${UBITRACK_COMPILE_DEFINITIONS} WINVER=0x0602)
  elseif( ${CMAKE_SYSTEM_VERSION} EQUAL 10.0 ) # Windows 10
    MESSAGE(STATUS "Setting minimum Windows version to Win8.1 WINVER=0x0603")
    set(UBITRACK_COMPILE_DEFINITIONS ${UBITRACK_COMPILE_DEFINITIONS} WINVER=0x0602)
  else() # Some other Windows
    MESSAGE(STATUS "Setting minimum Windows version to Vista WINVER=0x0600")
    set(UBITRACK_COMPILE_DEFINITIONS ${UBITRACK_COMPILE_DEFINITIONS} WINVER=0x0600)
  endif()
endif(WIN32)

IF(COMPILER_SUPPORTS_CXX11)
  set(UBITRACK_COMPILE_DEFINITIONS ${UBITRACK_COMPILE_DEFINITIONS} COMPILER_USE_CXX11)
ENDIF(COMPILER_SUPPORTS_CXX11)

IF(NOT WIN32)
  MESSAGE(STATUS "set boost::ublas alignment to 16")
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

