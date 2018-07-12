CPP = g++
FLAGS  = -march=native -O3 -funroll-loops


MemoryFill.o: MemoryFill.src/MemoryFill.cpp
	$(CPP) -c $(FLAGS) MemoryFill.src/MemoryFill.cpp -o MemoryFill.src/MemoryFill.o

StopWatch.o: CopMEM.src/StopWatch.cpp CopMEM.src/StopWatch.h
	$(CPP) -c $(FLAGS) CopMEM.src/StopWatch.cpp  -o CopMEM.src/StopWatch.o

CopMEM.o: CopMEM.src/CopMEM.cpp
	$(CPP) -c $(FLAGS) CopMEM.src/CopMEM.cpp  -o CopMEM.src/CopMEM.o		


city.o: CopMEM.src/city.cpp
	$(CPP) -c $(FLAGS) CopMEM.src/city.cpp  -o CopMEM.src/city.o		

xxhash.o: CopMEM.src/xxhash.c
	$(CPP) -c $(FLAGS) CopMEM.src/xxhash.c  -o CopMEM.src/xxhash.o

metrohash64.o:
	$(CPP) -c $(FLAGS) CopMEM.src/metrohash64.cpp  -o CopMEM.src/metrohash64.o


CopMEM: CopMEM.o StopWatch.o xxhash.o metrohash64.o city.o
	$(CPP) $(FLAGS) CopMEM.src/CopMEM.o CopMEM.src/StopWatch.o CopMEM.src/xxhash.o CopMEM.src/metrohash64.o CopMEM.src/city.o -o copmem

MemoryFill: MemoryFill.o city.o
	$(CPP) $(FLAGS) MemoryFill.src/MemoryFill.o -o memoryfill


clean:
	rm -R -f *.o copmem memoryfill

all: CopMEM MemoryFill
