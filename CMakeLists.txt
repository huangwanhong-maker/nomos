# nomos: the methods of the social sciences.
#
# Every module is a target nomos::<module>. nomos uses the targets of calx, morphe, physis and
# tyche: inside The Science Infrastructures, the siblings foundations/calx, foundations/morphe,
# foundations/physis and foundations/tyche (each added here when no larger project has added it);
# built alone elsewhere, they would be fetched at pinned versions (not yet written). Built on its
# own, nomos builds its checks (NOMOS_CHECKS); added to a larger project, only when asked.

cmake_minimum_required(VERSION 3.20)
project(nomos VERSION 0.1.0 LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(NOMOS_TOP_LEVEL ON)
else()
  set(NOMOS_TOP_LEVEL OFF)
endif()
option(NOMOS_CHECKS "Build nomos's checks" ${NOMOS_TOP_LEVEL})

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

# nomos_sibling(<directory> <a target it defines>): the sibling library, added once.
function(nomos_sibling directory target)
  if(NOT TARGET ${target})
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../${directory}/CMakeLists.txt)
      add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../${directory} ${CMAKE_BINARY_DIR}/${directory})
    else()
      message(FATAL_ERROR "nomos needs ${directory}: build it inside The Science Infrastructures (foundations/${directory}).")
    endif()
  endif()
endfunction()

# The siblings in their order: calx and morphe, then physis, then tyche.
nomos_sibling(calx calx::linalg)
nomos_sibling(morphe morphe::topology)
nomos_sibling(physis physis::units)
nomos_sibling(tyche tyche::probability)
if(NOMOS_CHECKS)
  enable_testing()
endif()

# The modules in their order (decision 0003, §2.2).
add_subdirectory(measurement)
add_subdirectory(indices)
add_subdirectory(games)
add_subdirectory(matching)
add_subdirectory(choice)
add_subdirectory(econometrics)
add_subdirectory(agents)
