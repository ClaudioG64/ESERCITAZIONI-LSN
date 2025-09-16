CXX = g++
CXXFLAGS = -Wall -O2

# target principale
all: 8_1.exe 8_2.exe

# link finale per 8_1
8_1.exe: 8_1.o random.o
	$(CXX) $(CXXFLAGS) -o $@ $^

8_1.o: 8_1.cpp random.h
	$(CXX) $(CXXFLAGS) -c 8_1.cpp

# link finale per 8_2
8_2.exe: 8_2.o random.o
	$(CXX) $(CXXFLAGS) -o $@ $^

8_2.o: 8_2.cpp random.h
	$(CXX) $(CXXFLAGS) -c 8_2.cpp

# random
random.o: random.cpp random.h
	$(CXX) $(CXXFLAGS) -c random.cpp

# pulizia
clean:
	rm -f *.o *.exe *.dat
