CC=g++ -std=c++11
#STUFF= -Wall -g -fopenmp -pthread -Ofast
STUFF= -Wall -g -O1
LDIR=-L../L2-lib/
LIB=-lboost_program_options -lboost_system -lboost_filesystem -lm -llpo
IDIR=-I../L2-inc/

OBJECTS=lomex2.o fun_kmers.o file_reader.o file_writer.o kmer_index.o fun_poa.o
TARGET=lomex2

clean:
	rm -f $(TARGET) $(OBJECTS)

lomex2: fun_kmers.o file_reader.o file_writer.o kmer_index.o fun_poa.o lomex2.o
	$(CC) -o $(TARGET) fun_kmers.o file_reader.o file_writer.o kmer_index.o fun_poa.o lomex2.o $(LDIR) $(LIB) $(STUFF)

poatest: main-poa.o
	$(CC) -o poatest main-poa.o $(LDIR) $(LIB) $(STUFF)

main-poa.o: main-poa.cpp
	$(CC) -c main-poa.cpp $(IDIR) $(STUFF)

fun_poa.o: fun_poa.cpp
	$(CC) -c fun_poa.cpp $(IDIR) $(STUFF)

lomex2.o: lomex2.cpp
	$(CC) -c lomex2.cpp $(STUFF) $(IDIR)

fun_kmers.o:
	$(CC) -c fun_kmers.cpp $(STUFF)


