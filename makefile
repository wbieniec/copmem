CPP = g++
FLAGS  = -march=native -O3 -funroll-loops


MemoryFill.o: MemoryFill.src/MemoryFill.cpp
	$(CPP) -c $(FLAGS) MemoryFill.src/MemoryFill.cpp -o MemoryFill.src/MemoryFill.o

StopWatch.o: CopMEM.src/StopWatch.cpp CopMEM.src/StopWatch.h
	$(CPP) -c $(FLAGS) CopMEM.src/StopWatch.cpp  -o CopMEM.src/StopWatch.o

CopMEM.o: CopMEM.src/CopMEM.cpp
	$(CPP) -c $(FLAGS) CopMEM.src/CopMEM.cpp  -o CopMEM.src/CopMEM.o		

CopMEM: CopMEM.o StopWatch.o
	$(CPP) $(FLAGS) CopMEM.src/CopMEM.o CopMEM.src/StopWatch.o -o copmem

MemoryFill: MemoryFill.o
	$(CPP) $(FLAGS) MemoryFill.src/MemoryFill.o -o memoryfill


clean:
	rm -R -f *.o copmem memoryfill

all: CopMEM MemoryFill
