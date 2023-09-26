
SRCDIR=./

CC = gcc

#  Not given `-DMKLROOT`.  Use own coded SVD/EVD functions .
## MKLROOT = /opt/intel/mkl
## MKLLIBDIR = $(MKLROOT)/lib/intel64
#CFLAGS = -g -Wall -Wextra -DDEBUG
CFLAGS = -g -Wall -Wextra
LIBS = -lm
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
## LIBS = -Wl,--start-group $(MKLLIBDIR)/libmkl_intel_lp64.a $(MKLLIBDIR)/libmkl_sequential.a $(MKLLIBDIR)/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

PROGS = msfit mmfit markpdb anisovar genconfms

all: $(PROGS)

msfit: msfit.o superimposems.o covarms.o pdb.o memory.o matrix.o chol.o svdgl.o lapack.o
	$(CC) -o $@ $^ -mcmodel=medium $(LDFLAGS) $(LIBS)



mmfit: mmfit.o superimpose.o mmvars.o pdb.o map.o memory.o matrix.o chol.o svdgl.o lapack.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)


anisovar: anisovar.o matrix.o chol.o
	$(CC) -o $@ $^ -mcmodel=medium $(LDFLAGS) $(LIBS)

markpdb: markpdb.o markpdb_utils.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)

markpdb.o: $(addprefix $(SRCDIR)/, markpdb.c molec.h atypes.h pdb.h pdbconst.h lapack.h)
	$(CC) -c -g -Wall -Wextra -DDEBUG -I$(SRCDIR) $< -o $@

genconfms: genconfms.o rand.o superimposems.o memory.o matrix.o markpdb_utils.o svdgl.o chol.o lapack.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)

clean::
	rm -f *.o

veryclean::
	rm -rf *.o ${PROGS}


%.o: $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) -I$(SRCDIR) $< -o $@

# markpdb.o: $(addprefix $(SRCDIR)/, molec.h atypes.h pdb.h pdbconst.h lapack.h)

lapack.o: $(addprefix $(SRCDIR)/, lapack.c atypes.h lapack.h)
	$(CC) $(CFLAGS) -Wno-unused-parameter -c $< -o $@
