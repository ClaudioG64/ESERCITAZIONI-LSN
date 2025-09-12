CXX = g++
CXXFLAGS = -Wall -O2

# Sorgenti degli esercizi
SRC = $(wildcard es*.cpp)

# Eseguibili (senza .cpp)
EXE = $(SRC:.cpp=)

# File comuni (sempre da compilare/linkare)
COMMON = random.cpp

all: $(EXE)

# Regola generica: compila ogni es1_X.cpp insieme a random.cpp
%: %.cpp $(COMMON)
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(EXE) *.o
