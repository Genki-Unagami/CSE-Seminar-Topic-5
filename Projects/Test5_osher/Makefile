# *********************************************************************************************
# README
# *********************************************************************************************
#
# Available configuration Parameters for ExaHyPE
#
# export variable  |  default-value  |  further values         |  description
#-------------------------------------------------------------------------------------------------------
# COMPILER            Intel             GNU                          Used compiler (and linker)
# MODE                Release           Debug,                       Verbosity and Debug level
#                                       Asserts,PeanoProfile,Profile          
# SHAREDMEM           None              OMP, TBB                     Shared-memory parallelisation
# DISTRIBUTEDMEM      None              MPI                          Distributed-memory parallelisation
# BOUNDARYCONDITIONS  None              Periodic                     Type of boundary conditions
# USE_IPO             Off               On, Unsafe                   IPO during compilation (intel only)
# ******************************************************************************************************

# Helper
SHELL = bash
tolower = $(shell echo $(1) | tr '[:upper:]' '[:lower:]')

# Default Paths and name
PEANO_KERNEL_PEANO_PATH=/home/genki/Desktop/CSE_seminar/ExaHyPE-Engine/./Peano/peano
PEANO_KERNEL_TARCH_PATH=/home/genki/Desktop/CSE_seminar/ExaHyPE-Engine/./Peano/tarch
PEANO_TOOLBOX_MULTISCALELINKEDCELL_PATH=/home/genki/Desktop/CSE_seminar/ExaHyPE-Engine/./Peano/multiscalelinkedcell
PEANO_TOOLBOX_SHAREDMEMORY_ORACLES_PATH=/home/genki/Desktop/CSE_seminar/ExaHyPE-Engine/./Peano/sharedmemoryoracles
PEANO_TOOLBOX_MPI_BLANCING_PATH=/home/genki/Desktop/CSE_seminar/ExaHyPE-Engine/./Peano/mpibalancing
EXAHYPE_PATH=/home/genki/Desktop/CSE_seminar/ExaHyPE-Engine/./ExaHyPE
PROJECT_PATH=/home/genki/Desktop/CSE_seminar/ExaHyPE-Engine/./Demonstrators/Test5_osher
EXECUTABLE=ExaHyPE-Euler

# Default flags from specification file constants
PROJECT_CFLAGS+= -DDim2
ARCHITECTURE=CPU

# Default sharedmem/distributedmem required by the specification file
ifeq ($(SHAREDMEM),)
    SHAREDMEM=None
endif
ifeq ($(DISTRIBUTEDMEM),)
    DISTRIBUTEDMEM=None
endif

# If using TBB only, link pthread (comment this if you want to link a specific pthread)
ifeq ($(call tolower,$(SHAREDMEM)),tbb)
    ifeq ($(call tolower,$(DISTRIBUTEDMEM)),none)
        PROJECT_LFLAGS+= -pthread
    endif
endif

# Default options required by the specification file



# Main Makefile include
-include /home/genki/Desktop/CSE_seminar/ExaHyPE-Engine/./ExaHyPE/Makefile


all: 
	@echo Euler
	@echo =================
	@echo An ExaHyPE solver