
cmake_minimum_required(VERSION 2.4)
#Change PROJECT_NAME to the name of your project
PROJECT(sccan)

OPTION(USE_ITK "Use ITK Libraries" ON)
set (CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
  
# Set up ITK
IF(USE_ITK)
  FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
      "Cannot build without ITK.  Please set ITK_DIR.")
  ENDIF(ITK_FOUND)
ENDIF(USE_ITK)

#The following lines are required to use Dart
INCLUDE(CTest)
ENABLE_TESTING()

INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake) 
IF (USE_ITK_FILE) 
  INCLUDE(${USE_ITK_FILE}) 
ENDIF(USE_ITK_FILE) 

SET(SCCAN_INCLUDE_DIRS
./
)
SET(SCCAN_SRC_DIR
./  
)

INCLUDE_DIRECTORIES(${SCCAN_INCLUDE_DIRS})
LINK_DIRECTORIES( ${ITK_LIBRARY_PATH} )

# non-templated class -- this should be stored in a library and linked in...
SET(SOURCES "./antsCommandLineParser" "./antsCommandLineOption" )

ADD_EXECUTABLE(sccan ${SCCAN_SRC_DIR}sccan.cxx ${SOURCES})
TARGET_LINK_LIBRARIES(sccan ITKCommon ITKIO ITKNumerics )
