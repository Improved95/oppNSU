cmake_minimum_required(VERSION 3.16.3)
project(cluster)
set(exe_name ${PROJECT_NAME})

if (MSVC)
    # warning level 4 and all warnings as errors
    add_compile_options(/W4 /WX)
else()
    # lots of warnings and all warnings as errors
    add_compile_options(-Wall -Wextra -pedantic -Werror)
endif()

find_package(MPI REQUIRED)
find_package(Threads REQUIRED)

add_library(task_queue STATIC task_queue.c)
add_executable(${exe_name} main.c)

include_directories(PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(${exe_name} PUBLIC task_queue PUBLIC ${MPI_C_LIBRARIES} PUBLIC Threads::Threads)
