CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2

TARGET = tsp_ga
SOURCES = main.cpp TSPGA.cpp random.cpp
HEADERS = City.h Individual.h TSPGA.h random.h

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES)

clean:
	rm -f $(TARGET)

.PHONY: clean