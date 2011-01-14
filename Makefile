CC=g++
CFLAGS=-Wall -pedantic -ggdb
OBJECTS=main.o RunParameter.o F1.o F2.o F3.o F4.o F5.o F6.o F7.o F8.o F9.o F10.o F11.o F12.o F13.o F14.o F15.o F16.o F17.o F18.o F19.o F20.o Benchmarks.o

main: $(OBJECTS)
	$(CC) $(CFLAGS) -o main $(OBJECTS)

main.o: Header.h RunParameter.h F1.h F2.h F3.h F4.h F5.h F6.h F7.h F8.h F9.h F10.h F11.h F12.h F13.h F14.h F15.h F16.h F17.h F18.h F19.h F20.h Benchmarks.h main.cpp
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

F5.o: F5.h Benchmarks.h F5.cpp
	$(CC) $(CFLAGS) -c F5.cpp

F6.o: F6.h Benchmarks.h F6.cpp
	$(CC) $(CFLAGS) -c F6.cpp

F7.o: F7.h Benchmarks.h F7.cpp
	$(CC) $(CFLAGS) -c F7.cpp

F8.o: F8.h Benchmarks.h F8.cpp
	$(CC) $(CFLAGS) -c F8.cpp

F9.o: F9.h Benchmarks.h F9.cpp
	$(CC) $(CFLAGS) -c F9.cpp

F10.o: F10.h Benchmarks.h F10.cpp
	$(CC) $(CFLAGS) -c F10.cpp

F11.o: F11.h Benchmarks.h F11.cpp
	$(CC) $(CFLAGS) -c F11.cpp

F12.o: F12.h Benchmarks.h F12.cpp
	$(CC) $(CFLAGS) -c F12.cpp

F13.o: F13.h Benchmarks.h F13.cpp
	$(CC) $(CFLAGS) -c F13.cpp

F14.o: F14.h Benchmarks.h F14.cpp
	$(CC) $(CFLAGS) -c F14.cpp

F15.o: F15.h Benchmarks.h F15.cpp
	$(CC) $(CFLAGS) -c F15.cpp
	
F16.o: F16.h Benchmarks.h F16.cpp
	$(CC) $(CFLAGS) -c F16.cpp

F17.o: F17.h Benchmarks.h F17.cpp
	$(CC) $(CFLAGS) -c F17.cpp
	
F18.o: F18.h Benchmarks.h F18.cpp
	$(CC) $(CFLAGS) -c F18.cpp
	
F19.o: F19.h Benchmarks.h F19.cpp
	$(CC) $(CFLAGS) -c F19.cpp

F20.o: F20.h Benchmarks.h F20.cpp
	$(CC) $(CFLAGS) -c F20.cpp

.PHONY : clean
clean:
	rm -f main output $(OBJECTS)
