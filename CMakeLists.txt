cmake_minimum_required(VERSION 3.3)
project(LWGA)

set(LWGA_CXX_FLAGS "-pthread -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic")
set(CMAKE_CXX_FLAGS "-std=c++14 ${LWGA_CXX_FLAGS}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ${LWGA_CXX_FLAGS}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -ggdb -std=c++14 ${LWGA_CXX_FLAGS}")

#set(SOURCE_FILES
#        include/cop_ga.h
#        src/cop_ga.cpp
#        src/cop_ga_main.cpp
#        src/ctop_ga.cpp
#    )

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

set(ICRA17_COP
        ${GENERIC_FILES}
        include/cop_ga.h
        src/cop_ga.cpp
        runs/icra17_cop.cpp)

set(PATH_RELINKING
        include/path_relinking.h)

set(TOP_GA
        ${GENERIC_FILES}
        include/top_ga.h
        src/top_ga_main.cpp
        src/top_ga.cpp
        src/ga_utils.cpp)

set(ICRA18_CTOP_MIQP_COMPARE
        ${GENERIC_FILES}
        include/ctop_ga.h
        src/ctop_ga.cpp
        runs/icra18_miqp_compare.cpp)

include_directories(.)
include_directories(include)

#add_executable(LWGA ${SOURCE_FILES})
add_executable(cop-ga ${COP_GA})
target_link_libraries(cop-ga pthread)

add_executable(ctop-ga ${CTOP_GA})
target_link_libraries(ctop-ga pthread)

add_executable(icra17-cop ${ICRA17_COP})
target_link_libraries(icra17-cop pthread)

add_executable(top_ga ${TOP_GA})
target_link_libraries(top_ga pthread)

add_executable(icra18_ctop_miqp_compare ${ICRA18_CTOP_MIQP_COMPARE})
target_link_libraries(icra18_ctop_miqp_compare pthread)