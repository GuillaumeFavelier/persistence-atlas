cmake_minimum_required(VERSION 3.2)
project(ttk VERSION 0.9.4 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

option(TTK_ENABLE_KAMIKAZE "Enable Kamikaze compilation mode" OFF)
option(TTK_ENABLE_CPU_OPTIMIZATION "Enable native CPU optimizations" ON)

set(TTK_BUILD_VTK_WRAPPERS ON)
set(TTK_BUILD_STANDALONE_APPS ON)

add_subdirectory(core)
add_subdirectory(standalone)

ttk_print_summary()
