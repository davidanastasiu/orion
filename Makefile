######
# Makefile for Orion+ by David C. Anastasiu
# Dependencies: libm, libsvd, libl2ap, libGKlib
######

###
# Compile choices
###

# Define libraries and library directories
OMP := -fopenmp
LIBDIRS := -L"${HOME}/programs/lib" -L/usr/lib -L/usr/local/lib
LIBS := -lm -lsvd -lGKlib -ll2ap
# Source include directories
INC += -I"${HOME}/programs/include" -I/usr/local/include/
# C flags
#ifdef DEBUG
    OPTOPTS := -g 
#else
#    OPTOPTS := -march=native -O3 -DNDEBUG
#endif
CFLAGS += -c -fno-strict-aliasing -DLINUX -D_FILE_OFFSET_BITS=64 -std=c99 -fPIC -Wall -Wno-unused-function -Wno-unused-label -Wno-unused-variable -Wno-parentheses -Wsequence-point -Wno-unknown-pragmas
# Other compile choices
EXENAME := orion+
BDIR := build
SDIR := src
RM := rm -rf
CC := gcc
CP := cp

###
# Execute make
###

# Add inputs and outputs from these tool invocations to the build variables 
HEADERS := $(shell ls $(SDIR)/*.h)
C_SRCS := $(shell cd $(SDIR) && ls *.c)
C_OBJS := $(C_SRCS:%.c=$(BDIR)/%.o)

# All Targets
all: $(EXENAME) 


# Program depends on objects
$(EXENAME): $(C_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	$(CC) $(LIBDIRS) $(OMP) -o $(EXENAME) $(C_OBJS) $(OMP) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '
	
# Objects depend on its source and all headers
$(BDIR)/%.o: $(SDIR)/%.c $(HEADERS)
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	$(CC) $(INC) $(OMP) $(OPTOPTS) $(CFLAGS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# Clean Target
clean:
	$(RM) $(BDIR)/*.o $(BDIR)/*.d $(BDIR)/$(EXENAME) $(EXENAME)
	@echo ' '

# These targets do not produce files
.PHONY: clean install

