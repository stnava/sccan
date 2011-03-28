project(sccan)
cmake_minimum_required(VERSION 2.6.2)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} )
find_package(Eigen3 REQUIRED)
include_directories(EIGEN3_INCLUDE_DIR)

# set (CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

 
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


INCLUDE_DIRECTORIES(${SUMBA_INCLUDE_DIRS})


# /Users/brianavants/code/eigen/
set(CMAKE_BUILD_TYPE "Release")

SET(SCCAN_INCLUDE_DIRS
${EIGEN3_INCLUDE_DIR}
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
TARGET_LINK_LIBRARIES(sccan ${ITK_LIBRARIES} )
