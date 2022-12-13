# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(BLoM)

if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(blom blom.C)
target_link_libraries(blom Garfield::Garfield)

# ---Build executable------------------------------------------------------------
add_executable(gen_gas gen_gas.cpp)
target_link_libraries(gen_gas Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
foreach(_file data/BLoM/BLoM_dielectrics.dat data/BLoM/BLoM_field.txt data/BLoM/BLoM_ramo.txt data/BLoM/BLoM.mphtxt)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()
