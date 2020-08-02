cmake_minimum_required(VERSION 3.3)
project(LWGA)

set(CMAKE_CXX_STANDARD 14)
set(LWGA_CXX_FLAGS "-pthread -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wconversion")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LWGA_CXX_FLAGS}")
#set(CMAKE_CXX_FLAGS "-std=c++14 ${LWGA_CXX_FLAGS}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ${LWGA_CXX_FLAGS}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -ggdb -std=c++14 ${LWGA_CXX_FLAGS}")

#set(SOURCE_FILES
#        include/cop_ga.h
#        src/cop_ga.cpp
#        src/cop_ga_main.cpp
#        src/ctop_ga.cpp
#    )

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

set(GENERIC_FILES
        include/ga_types.h include/ga_utils.h src/ga_utils.cpp)

set(COP_GA
        ${GENERIC_FILES}
        include/cop_ga.h
        src/cop_ga.cpp
        src/cop_ga_main.cpp)

set(CTOP_GA
        ${GENERIC_FILES}
        include/ctop_ga.h
        src/ctop_ga.cpp
        src/ctop_ga_main.cpp)

set(CTOP_TUNING
        ${GENERIC_FILES}
        include/ctop_ga.h
        include/tuning_instances.h
        src/ctop_ga.cpp
        src/crs_tuner.cpp)

set(ICRA17_COP
        ${GENERIC_FILES}
        include/cop_ga.h
        src/cop_ga.cpp
        runs/icra17_cop.cpp)

set(PATH_RELINKING
        include/path_relinking.h)

#set(TOP_GA
#        ${GENERIC_FILES}
#        include/top_ga.h
#        src/top_ga_main.cpp
#        src/top_ga.cpp
#        src/ga_utils.cpp)

set(ICRA18_CTOP_MIQP_COMPARE
        ${GENERIC_FILES}
        include/ctop_ga.h
        src/ctop_ga.cpp
        runs/icra18_miqp_compare.cpp)

include_directories(.)
include_directories(include)
include_directories(/opt/gurobi811/linux64/include)

link_directories(/opt/gurobi811/linux64/lib)

add_library(dubins-curves src/dubins.cpp)

add_library(ga-utils src/ga_utils.cpp)
target_link_libraries(ga-utils dubins-curves)

add_library(rng src/rng.cpp)

add_library(mdmtsp src/mdmtsp_solver.cpp)
target_link_libraries(mdmtsp gurobi_c++ gurobi81 pthread)

add_library(cop_ga src/cop_ga.cpp)
target_link_libraries(cop_ga ga-utils pthread)

add_library(ctop_ga src/ctop_ga.cpp)
target_link_libraries(ctop_ga ga-utils pthread)

add_library(two_level_solver src/two_level_solver.cpp)
target_link_libraries(two_level_solver mdmtsp cop_ga)

add_library(pop_ga src/pop_ga.cpp)
target_link_libraries(pop_ga ga-utils rng)

add_library(libtsp_reader src/libtsp_reader.cpp)

#add_executable(LWGA ${SOURCE_FILES})
add_executable(cop-ga ${COP_GA})
target_link_libraries(cop-ga dubins-curves pthread)

add_executable(ctop-ga ${CTOP_GA})
target_link_libraries(ctop-ga dubins-curves pthread)

add_executable(icra17-cop ${ICRA17_COP})
target_link_libraries(icra17-cop dubins-curves pthread)

add_executable(pop-ga src/pop_ga_main.cpp)
target_link_libraries(pop-ga pop_ga libtsp_reader pthread)

#add_executable(top_ga ${TOP_GA})
#target_link_libraries(top_ga pthread)

add_executable(icra18_ctop_miqp_compare ${ICRA18_CTOP_MIQP_COMPARE})
target_link_libraries(icra18_ctop_miqp_compare dubins-curves pthread)

add_executable(crs_tuner ${CTOP_TUNING})
target_link_libraries(crs_tuner dubins-curves pthread)

add_library(dcop-ga src/dubins.cpp src/ga_utils.cpp src/dcop_ga.cpp)
target_link_libraries(dcop-ga pthread)

add_executable(dubins-crs-tuner src/dcop_crs_tuner.cpp)
target_link_libraries(dubins-crs-tuner dcop-ga)

add_executable(dcop-ga-main  src/dcop_ga_main.cpp)
target_link_libraries(dcop-ga-main dcop-ga)

add_executable(iros18-experiments-rho runs/iros18_experiments_rho.cpp)
target_link_libraries(iros18-experiments-rho dcop-ga)

add_executable(iros18-experiments-cc runs/iros18_experiments_cc.cpp)
target_link_libraries(iros18-experiments-cc dcop-ga)

add_executable(mdmtsp_problem src/mdmtsp_main.cpp src/ga_utils.cpp)
target_link_libraries(mdmtsp_problem dubins-curves mdmtsp)

add_executable(two-level-solver-main src/two_level_solver_main.cpp)
target_link_libraries(two-level-solver-main two_level_solver)

add_executable(large-ctop src/large_ctop_main.cpp)
target_link_libraries(large-ctop ctop_ga)

add_executable(cop-crs-tuner src/cop_crs_tuner.cpp)
target_link_libraries(cop-crs-tuner cop_ga)

add_executable(multi-area-comparisson src/multi-area-comparisson.cpp)
target_link_libraries(multi-area-comparisson two_level_solver ctop_ga)

add_executable(chap7 runs/chap-7-runs.cpp)
target_link_libraries(chap7 two_level_solver ctop_ga)

add_executable(dop_runs runs/dop_runs.cpp)
target_link_libraries(dop_runs dcop-ga)

add_executable(dcop_online_replan runs/dcop_online_replan.cpp)
target_link_libraries(dcop_online_replan dcop-ga)

include(CTest)
if (BUILD_TESTING)
  add_subdirectory(test)
endif()
