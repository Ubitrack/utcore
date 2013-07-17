# ----------------------------------------------------------------------------
#  CMake file for custom boost bindings. See root CMakeLists.txt
#
# ----------------------------------------------------------------------------

# All header files for BOOSTBINDINGS
FILE(GLOB BOOSTBINDINGS_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "${CMAKE_CURRENT_SOURCE_DIR}/boost/numeric/bindings/*/*.h"  "${CMAKE_CURRENT_SOURCE_DIR}/boost/numeric/bindings/*/*.h"  
   "${CMAKE_CURRENT_SOURCE_DIR}/boost/numeric/bindings/*/*.hpp"  "${CMAKE_CURRENT_SOURCE_DIR}/boost/numeric/bindings/*/*.hpp"  
   "${CMAKE_CURRENT_SOURCE_DIR}/boost/numeric/bindings/*/*/*.hpp"  "${CMAKE_CURRENT_SOURCE_DIR}/boost/numeric/bindings/*/*/*.hpp"  
)

foreach(hdr ${BOOSTBINDINGS_HEADERS})
   string(REGEX REPLACE "^.*modules/utcore/3rd/boost-bindings/" "" hdr2 "${hdr}")
   GET_FILENAME_COMPONENT(fpath ${hdr2} PATH)
   #MESSAGE(STATUS "${UBITRACK_INCLUDE_INSTALL_PATH}/${fpath}")
   IF(fpath)
     install(FILES ${hdr} DESTINATION "${UBITRACK_INCLUDE_INSTALL_PATH}/${fpath}" COMPONENT main)
   ELSE(fpath)
     install(FILES ${hdr} DESTINATION "${UBITRACK_INCLUDE_INSTALL_PATH}" COMPONENT main)
   ENDIF(fpath)
endforeach()

