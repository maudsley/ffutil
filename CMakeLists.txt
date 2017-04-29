cmake_minimum_required(VERSION 3.0)

# where is boost?
include_directories(/opt/boost)
link_directories(/opt/boost/stage/lib)

# where is gtest?
include_directories(/opt/googletest/googletest/include)
link_directories(/opt/googletest/googlemock/gtest)

# enable C++11 and warnings
set(CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "-Wall ${CMAKE_CXX_FLAGS}")

include_directories(.)

# build unit test executable

set(TestSources
  finite_field.cpp
  unit_test.cpp
)

set(TestLibs
  gtest
  pthread
)

add_executable(unit_test ${TestSources})
target_link_libraries(unit_test ${TestLibs})
