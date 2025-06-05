.PHONY: all clean sample1 sample2 sample3 sample4

PREFIX=..
CC=mpiicc
FC=mpiifort
CFLAGS=-I${PREFIX}/include
LDFLAGS=-L${PREFIX}/lib
LIBS=-lctca
F_LIBS=-lctca_f -lctca_mod

SAMPLE1_PROGS=sample1_requester sample1_coupler sample1_worker1 sample1_worker2 sample1_worker3
SAMPLE2_PROGS=sample2_requester sample2_coupler sample2_worker1 sample2_worker2 sample2_worker3
SAMPLE3_PROGS=sample3_requester sample3_coupler sample3_worker
SAMPLE4_PROGS=sample4_requester sample4_coupler

all: sample1 sample2 sample3 sample4

sample1: $(SAMPLE1_PROGS)
sample2: $(SAMPLE2_PROGS)
sample3: $(SAMPLE3_PROGS)
sample4: $(SAMPLE4_PROGS)

%: %.c
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@
%: %.f90
	$(FC) $(CFLAGS) $(LDFLAGS) $^ $(F_LIBS) $(LIBS) -o $@

clean:
	/bin/rm -f *.o *.mod $(SAMPLE1_PROGS) $(SAMPLE2_PROGS) $(SAMPLE3_PROGS) $(SAMPLE4_PROGS)
