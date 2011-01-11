CC=g++
CFLAGS=-Wall -pedantic -ggdb
OBJECTS=main.o RunParameter.o F1.o F2.o F3.o F4.o Benchmarks.o

main: $(OBJECTS)
	$(CC) $(CFLAGS) -o main $(OBJECTS)

main.o: Header.h RunParameter.h F1.h F2.h F3.h F4.h Benchmarks.h main.cpp
	$(CC) $(CFLAGS) -c main.cpp 

Benchmarks.o: RunParameter.h Benchmarks.cpp
	$(CC) $(CFLAGS) -c Benchmarks.cpp

RunParameter.o: RunParameter.h Header.h RunParameter.cpp
	$(CC) $(CFLAGS) -c RunParameter.cpp

F1.o: F1.h Benchmarks.h F1.cpp
	$(CC) $(CFLAGS) -c F1.cpp

F2.o: F2.h Benchmarks.h F2.cpp
	$(CC) $(CFLAGS) -c F2.cpp
	
F3.o: F3.h Benchmarks.h F3.cpp
	$(CC) $(CFLAGS) -c F3.cpp

F4.o: F4.h Benchmarks.h F4.cpp
	$(CC) $(CFLAGS) -c F4.cpp

.PHONY : clean
clean:
	rm -f main $(OBJECTS)
