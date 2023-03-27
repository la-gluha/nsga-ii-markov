cmake_minimum_required(VERSION 3.19)
project(nsga_II__c_code C)

set(CMAKE_C_STANDARD 11)

include_directories(nsga2code)

add_executable(nsga_II__c_code
        nsga2code/crossover.h
        nsga2code/decode.h
        nsga2code/func-con.h
        nsga2code/init.h
        nsga2code/input.h
        nsga2code/keepaliven.h
        nsga2code/mut.h
        nsga2code/nsga2.c
        nsga2code/rancon.h
        nsga2code/random.h
        nsga2code/ranking.h
        nsga2code/realcross2.h
        nsga2code/realinit.h
        nsga2code/realmut1.h
        nsga2code/report.h
        nsga2code/roullette.h
        nsga2code/select.h
        nsga2code/uniformxr.h nsga2code/mexTest.cpp nsga2code/helloMex.cpp nsga2code/nsga2_matlab.cpp)
