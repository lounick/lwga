cmake_minimum_required(VERSION 3.3)
project(LWGA)

set(LWGA_CXX_FLAGS " -pthread -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 ${LWGA_CXX_FLAGS}")

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

include_directories(.)
include_directories(include)

#add_executable(LWGA ${SOURCE_FILES})
add_executable(cop-ga ${COP_GA})
target_link_libraries(cop-ga pthread)

add_executable(ctop-ga ${CTOP_GA})
target_link_libraries(ctop-ga pthread)