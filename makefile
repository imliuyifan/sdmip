# Makefile definitions for CPLEX 12.4
#    16 May 2013
# For Ubuntu 12.04 LTS
# Written by Yifan
#------------------------------------------------------------


#------------------------------------------------------------
# System Description for My Computer
#------------------------------------------------------------
SYSTEM     = x86-64_sles10_4.1
LIBFORMAT  = static_pic

#------------------------------------------------------------
# Location of CPLEX installations
#------------------------------------------------------------
CPLEXDIR = /opt/ibm/ILOG/CPLEX_Studio124/cplex

#------------------------------------------------------------
# Compiler Selection (CC is for c)
#------------------------------------------------------------
CC  = gcc -O3
# CC = gcc -g
# ---------------------------------------------------------------------
# Compiler options
# ---------------------------------------------------------------------
COPT  = -m64 -fPIC -fno-strict-aliasing 

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------
CLNFLAGS  = -lcplex -lm -pthread

# CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
# Let's automate this process
CPLEXLIBDIR := $(shell find /opt ~ -name libcplex.a | cat > ~/cplexPath && cat ~/cplexPath | grep -m1 -E 'cplex/lib' | awk -F "/libcplex.a" '{print $$1}')

# CPLEXINCDIR   = $(CPLEXDIR)/include/ilcplex
CPLEXINCDIR := $(shell find /opt ~ -name cplex.h | cat > ~/cplexPath && cat ~/cplexPath | grep -m1 -E 'cplex/include' | awk -F "/cplex.h" '{print $$1}')

# ---------------------------------------------------------------------
# To execute the CPLEX examples.
# ---------------------------------------------------------------------

CFLAGS  = $(COPT) -I$(CPLEXINCDIR) -L$(CPLEXLIBDIR)

#------------------------------------------------------------

EXECUTABLE = sd
SOURCES = sd.c cuts.c sigma.c delta.c omega.c lambda.c utility.c theta.c \
prob.c cell.c soln.c improve.c solverc.c master.c subprob.c \
rvgen.c input.c parser.c supomega.c optimal.c log.c \
testout.c memory.c quad.c batch.c
HEADERS =cuts.h sigma.h delta.h omega.h lambda.h utility.h theta.h \
prob.h cell.h soln.h improve.h solver.h master.h subprob.h \
rvgen.h input.h parser.h supomega.h optimal.h log.h \
testout.h memory.h quad.h batch.h sdconstants.h sdglobal.h
OBJECTS = sd.o cuts.o sigma.o delta.o omega.o lambda.o utility.o theta.o \
prob.o cell.o soln.o improve.o solverc.o master.o subprob.o \
rvgen.o input.o parser.o supomega.o optimal.o log.o \
testout.o memory.o quad.o batch.o


# ------------------------------------------------------------

# make definitions:
%.o : %.c
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXECUTABLE) $(CFLAGS) $(CLNFLAGS)

$(OBJECTS) : $(HEADERS)
	$(CC) -c $(SOURCES) $(CFLAGS)

clean :
	rm -rf sd *.o *.out *.txt *~ # *.class
	rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp *.dat Summary *.sol
	rm -rf Batch_dual Master_Dual Batch_Obj Master_Obj Batch_x incumb
