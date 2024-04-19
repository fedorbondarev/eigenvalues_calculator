cmake_minimum_required(VERSION 3.25)
project(eigenvalues_calculator C)

set(CMAKE_C_STANDARD 99)

add_executable(eigenvalues_calculator
        main.c
        return_codes.h)
