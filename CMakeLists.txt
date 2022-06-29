cmake_minimum_required(VERSION 3.5)

project(mtgso12 LANGUAGES CXX Fortran)

if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/gama/)
    message(FATAL_ERROR
        "Download projet GNU Gama : "
        "git clone https://git.savannah.gnu.org/git/gama.git")
endif()

enable_language(Fortran)            # case sensitive name

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(gama/lib .)

add_executable(cmp-gso cmp-gso.cpp gso/mtgso1.f gso/mtgso2.f
               gama/lib/gnu_gama/adj/icgs.h gama/lib/gnu_gama/adj/icgs.cpp)

add_executable(lauchli lauchli.cpp gama/lib/gnu_gama/adj/icgs.h gama/lib/gnu_gama/adj/icgs.cpp gso/mtgso1.f)

add_executable(randmat randmat.cpp gama/lib/gnu_gama/adj/icgs.h gama/lib/gnu_gama/adj/icgs.cpp gso/mtgso1.f)
