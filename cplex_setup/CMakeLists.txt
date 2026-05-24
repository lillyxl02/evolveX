cmake_minimum_required(VERSION 3.20)
project(TSP_Optimization VERSION 1.0)

include_directories(/Applications/CPLEX_Studio2211/cplex/include)
include_directories(/Applications/CPLEX_Studio2211/concert/include)
include_directories(include)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

set(CPLEXDIR /Applications/CPLEX_Studio2211/cplex)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I. -I${CPLEXDIR}/include/ilcplex")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${CPLEXDIR}/lib/arm64_osx/static_pic -lcplex -lpthread -ldl")

add_executable(TSP_Optimization
    cplex_setup/src/main.c
    cplex_setup/src/algoSelector.c
    cplex_setup/src/parser.c
    cplex_setup/src/tsp.c
    cplex_setup/src/menu.c

)

target_link_libraries(TSP_Optimization PRIVATE dl m)
include(CPack)


