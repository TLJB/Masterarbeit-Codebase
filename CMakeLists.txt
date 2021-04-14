CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

FIND_PACKAGE(deal.II 9.0.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

SET(TARGET "main")
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# ADD_SUBDIRECTORY("fem")
INCLUDE_DIRECTORIES(
    "fem/include"
    "material/include"
    )
FILE(GLOB_RECURSE FEM_SOURCES fem/source/*.cpp)
SET(TARGET_SRC ${FEM_SOURCES})

PROJECT(${TARGET})

DEAL_II_INVOKE_AUTOPILOT()