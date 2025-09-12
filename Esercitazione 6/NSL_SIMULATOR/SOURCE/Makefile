CC = g++
# Cambia da C++11 a C++14
CFLAGS = -O3 --std=c++14
INCLUDE = -I/opt/homebrew/include
LIBRARY = -L/opt/homebrew/lib
AFLAGS = -larmadillo

all: simulator.exe

simulator.exe : system.o NSL_SIMULATOR.o particle.o random.o
	$(CC) system.o particle.o random.o NSL_SIMULATOR.o -o simulator.exe $(LIBRARY) $(AFLAGS)

NSL_SIMULATOR.o : NSL_SIMULATOR.cpp system.h particle.h
	$(CC) -c NSL_SIMULATOR.cpp -o NSL_SIMULATOR.o $(CFLAGS) $(INCLUDE)

system.o : system.cpp system.h random.h particle.h
	$(CC) -c system.cpp -o system.o $(CFLAGS) $(INCLUDE)

particle.o : particle.cpp particle.h
	$(CC) -c particle.cpp -o particle.o $(CFLAGS) $(INCLUDE)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS) $(INCLUDE)

clean :
	rm *.o simulator.exe

remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*