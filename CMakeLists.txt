cmake_minimum_required(VERSION 3.20)
project(pvmls_template_generator VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(CGAL REQUIRED)
find_package(Boost REQUIRED COMPONENTS graph)
find_package(JPEG REQUIRED)

add_library(pvmls_generator src/generator.cpp)
target_include_directories(pvmls_generator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(pvmls_generator PUBLIC CGAL::CGAL Boost::graph JPEG::JPEG)
target_compile_features(pvmls_generator PUBLIC cxx_std_20)

if(MSVC)
  target_compile_options(pvmls_generator PRIVATE /W4 /permissive-)
else()
  target_compile_options(pvmls_generator PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
endif()

add_executable(pvmls-template-generator src/main.cpp)
target_link_libraries(pvmls-template-generator PRIVATE pvmls_generator)

include(CTest)
if(BUILD_TESTING)
  add_executable(pvmls-tests tests/test_generator.cpp)
  target_link_libraries(pvmls-tests PRIVATE pvmls_generator JPEG::JPEG)
  add_test(NAME pvmls-generator-tests COMMAND pvmls-tests)
endif()
