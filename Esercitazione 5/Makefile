CC = g++
CFLAGS = -Wall -O3 --std=c++11
INCLUDES = -I/opt/homebrew/include
LIBS = -L/opt/homebrew/lib -larmadillo

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe $(LIBS)
	
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS) $(INCLUDES)
	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS) $(INCLUDES)

clean :
	rm -f *.o main.exe seed.out
