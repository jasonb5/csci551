CC := gcc

MPI_CC := mpicc

CPPFLAGS := -g -Wall -std=c99

LDLIBS := -lm -lmpi

all: analysis num_integration

analysis: analysis.c
	$(CC) $(CPPFLAGS) -lm -o $@ $<	

num_integration: num_integration.c
	$(MPI_CC) $(CPPFLAGS) $(LDLIBS) -o $@ $<

clean:
	rm -rf analysis
	rm -rf num_integration
