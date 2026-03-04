find_package(deal.II 9.5.0 REQUIRED  
  HINTS  
    ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR} ~/bin/deal.II  
)  
set(TARGET "solver")  
set(TARGET_SRC  
  ${TARGET}.cc  
)  
cmake_minimum_required(VERSION 3.13.4)  
deal_ii_initialize_cached_variables()  
project(${TARGET} CXX)  
deal_ii_invoke_autopilot()
