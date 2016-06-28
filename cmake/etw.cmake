# Event Tracing for Windows setup
#


# Check if OS supports ETW
MACRO(CHECK_ETW)
#  CMAKE_MC_COMPILER - where to find mc.exe
if (WIN32)
  # cmake has CMAKE_RC_COMPILER, but no message compiler
  if ("${CMAKE_GENERATOR}" MATCHES "Visual Studio")
    # this path is only present for 2008+, but we currently require PATH to
    # be set up anyway
    get_filename_component(sdk_dir "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Microsoft SDKs\\Windows;CurrentInstallFolder]" REALPATH)
    get_filename_component(kit_dir "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows Kits\\Installed Roots;KitsRoot]" REALPATH)
    if (X64)
      set(sdk_bindir "${sdk_dir}/bin/x64")
      set(kit_bindir "${kit_dir}/bin/x64")
    else (X64)
      set(sdk_bindir "${sdk_dir}/bin")
      set(kit_bindir "${kit_dir}/bin/x86")
    endif (X64)
  endif ()
  find_program(CMAKE_MC_COMPILER mc.exe HINTS "${sdk_bindir}" "${kit_bindir}"
    DOC "path to message compiler")

  # This should probably check for the Windows SDK as well here ..
    
  if (NOT CMAKE_MC_COMPILER)
     SET(ENABLE_ETW OFF CACHE BOOL "Event Tracing for Windows disabled.")
   ELSE()
     SET(ENABLE_ETW ON CACHE BOOL "Enable Event Tracing for Windows")
  endif (NOT CMAKE_MC_COMPILER)
  message(STATUS "Found message compiler: ${CMAKE_MC_COMPILER}")
  mark_as_advanced(CMAKE_MC_COMPILER)
endif(WIN32)

 SET(HAVE_ETW ${ENABLE_ETW})
ENDMACRO()

CHECK_ETW()

# Create provider headers
IF(ENABLE_ETW)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/tracing/etw/ubitrack_etw_providers.man.base
    ${CMAKE_BINARY_DIR}/utUtil/ubitrack_etw_providers.man COPYONLY)

 ADD_CUSTOM_COMMAND(
   OUTPUT "${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_etw.h" 
           "${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_etw.rc" 
           "${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_etwTemp.bin"
   COMMAND ${CMAKE_MC_COMPILER} -um ${CMAKE_BINARY_DIR}/utUtil/ubitrack_etw_providers.man -z probes_ubitrack_etw
   WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/utUtil"
   COMMENT "Exectuing: mc.exe -um ${CMAKE_BINARY_DIR}/utUtil/ubitrack_etw_providers.man -z probes_ubitrack_etw"
   DEPENDS ${CMAKE_BINARY_DIR}/utUtil/ubitrack_etw_providers.man
 )
  
ADD_CUSTOM_TARGET(gen_etw_header
  DEPENDS  
  ${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_etw.h
  ${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_etw.rc
  ${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_etwTemp.bin
  )

  install(FILES ${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_etw.h DESTINATION "${UBITRACK_INCLUDE_INSTALL_PATH}/utUtil" COMPONENT main)
ENDIF()

FUNCTION(ETW_INSTRUMENT target)
  ADD_DEPENDENCIES(${target} gen_etw_header)
ENDFUNCTION()