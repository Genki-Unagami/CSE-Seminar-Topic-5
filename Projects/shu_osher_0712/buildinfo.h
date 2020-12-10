/* This file is autogenerated by generate-buildinfo.sh */

#ifndef _EXAHYPE_BUILD_INFO_H_
#define _EXAHYPE_BUILD_INFO_H_

#define EXAHYPE_BUILDINFO_AVAILABLE
#define EXAHYPE_BUILD_DATE          "Mon 07 Dec 2020 10:04:23 AM CET"
#define EXAHYPE_BUILD_HOST          "genki-MacBookPro"

/* Strings passed by the Makefile */
#define EXAHYPE_BUILD_INFO \
    "COMPILER = GNU\n" \
    "MODE = Release\n" \
    "SHAREDMEM = None\n" \
    "DISTRIBUTEDMEM = None\n" \
    "------\n" \
    "ARCHITECTURE = CPU\n" \
    "CC = g++\n" \
    "BOUNDARYCONDITIONS = None\n" \
    "------\n" \
    "COMPILER_CFLAGS = -DTrackGridStatistics -std=c++11 -pedantic -Wall -Drestrict=__restrict__ -pipe -D__assume_aligned=__builtin_assume_aligned -Wstrict-aliasing -fopenmp-simd -O3 \n" \
    "COMPILER_LFLAGS = -lm -lstdc++ -lrt\n" \
    "FCOMPILER_CFLAGS = -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -cpp -Wall -Wno-tabs -Wno-unused-variable -fopenmp-simd -O3 \n" \
    "FCOMPILER_LFLAGS = \n" \
    "------\n" \
    "EXAHYPE_PATH = /home/genki/Desktop/CSE_seminar/ExaHyPE-Engine/./ExaHyPE\n" \
    "PROJECT_PATH = /home/genki/Desktop/CSE_seminar/ExaHyPE-Engine/./Demonstrators/entropy_osher_0712\n" \
    "PEANO_KERNEL_PEANO_PATH = /home/genki/Desktop/CSE_seminar/ExaHyPE-Engine/./Peano/peano\n" \
    "PEANO_KERNEL_TARCH_PATH = /home/genki/Desktop/CSE_seminar/ExaHyPE-Engine/./Peano/tarch\n" \
    "------\n" \
    "PROJECT_CFLAGS = -DDim2\n" \
    "PROJECT_LFLAGS = \n" \
    ""


/*
 * ExaHyPE Git Repository information extraction
 * Extracted from git repository in /home/genki/Desktop/CSE_seminar/ExaHyPE-Engine/./ExaHyPE
 */

/* Information collected with git version 2.25.1 */
#define EXAHYPE_GIT_INFO "master  31b5eb94b Fri Nov 6 13:10:00 2020"

/*
 * Peano Git Repository information extraction
 * Extracted from git repository in /home/genki/Desktop/CSE_seminar/ExaHyPE-Engine/./Peano/peano
 */

/* Information collected with git version 2.25.1 */
#define PEANO_GIT_INFO "HEAD  40b5c080 Tue Oct 13 22:36:28 2020"
#define PEANO_SVN_INFO  PEANO_GIT_INFO /* transition time */

/*
    Peano version check

    FIXME TODO This is the worst place ever to hook in version
               requirements. Please move the following to somewhere
               suitable, such as the Peano startup phase.

*/

#include "peano/version.h"

#endif /* EXAHYPE_BUILD_INFO_H */
#if PEANO_VERSION<2509 
#error Old Peano version. Version 2509 required. Please update your Peano installation.
#endif
