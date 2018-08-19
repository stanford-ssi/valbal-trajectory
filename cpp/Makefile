TARGET = traj
LIBS = -lm
CC = g++
CXXFLAGS = -g -Werror -Ofast -fno-strict-aliasing -Wall

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.h)

%.o: %.cpp $(HEADERS)
	$(CC) $(CXXFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@

diff: diff.o
	$(CC) diff.o -o difftest

clean:
	-rm -f *.o
	-rm -f $(TARGET)

open:
	vi -p traj.cpp config.h data.h data.cpp sim.h sim.cpp diff.h diff.cpp

opendiff:
	vi -p diff.h diff.cpp
