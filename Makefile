CC=mpicc

totient: TotientRange.c
	$(CC) -o totient TotientRange.c