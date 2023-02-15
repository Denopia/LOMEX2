CC=g++ -std=c++11
#STUFF= -Wall -g -fopenmp -pthread -Ofast
STUFF= -Wall -g -O3
#LDIR=-L../LII-lib/
LDIR=-L../spoa/build/lib/
#LIB=-lboost_program_options -lboost_system -lboost_filesystem -lm -llpo -fopenmp
#LIB=-lboost_program_options -lboost_system -lboost_filesystem -lm -lspoa -fopenmp
LIB= -lm -lspoa -fopenmp
#IDIR=-I../LII-inc/
IDIR=-I../spoa/include/spoa/
OMP=-fopenmp
OBJECTS=sake.o fun_kmers.o file_reader.o file_writer.o kmer_index.o fun_spoa.o
TARGET=sake

clean:
	rm -f $(TARGET) $(OBJECTS)

sake: fun_kmers.o file_reader.o file_writer.o kmer_index.o fun_spoa.o sake.o
	$(CC) -o $(TARGET) fun_kmers.o file_reader.o file_writer.o kmer_index.o fun_spoa.o sake.o $(LDIR) $(LIB) $(STUFF)

poatest: main-poa.o
	$(CC) -o poatest main-poa.o $(LDIR) $(LIB) $(STUFF)

main-poa.o: main-poa.cpp
	$(CC) -c main-poa.cpp $(IDIR) $(STUFF)

fun_spoa.o: fun_spoa.cpp
	$(CC) -c fun_spoa.cpp $(IDIR) $(STUFF)

sake.o: sake.cpp
	$(CC) -c sake.cpp $(STUFF) $(IDIR) $(OMP)

fun_kmers.o:
	$(CC) -c fun_kmers.cpp $(STUFF)


