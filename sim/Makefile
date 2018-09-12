TARGET = traj
LIBS = -lm -ladept -lpthread
CC = g++
CXXFLAGS = -g -Werror -Ofast -march=native -fno-strict-aliasing -Wall -ffast-math

PATHS = -I ../balloons-VALBAL/lib/spa
PATHS += -I ../balloons-VALBAL/src
PATHS += -I /usr/local/include


.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.h)

#Jank additions to compile lasagna and Pastasim. Must run utils/flight_code.sh first
OBJECTS += ../balloons-VALBAL/src/LasagnaController.cpp
HEADERS += ../balloons-VALBAL/src/LasagnaController.h
OBJECTS += ../balloons-VALBAL/src/Utils.cpp
HEADERS += ../balloons-VALBAL/src/Utils.h
OBJECTS += ../balloons-VALBAL/lib/spa/spa.cpp
HEADERS += ../balloons-VALBAL/lib/spa/spa.h
OBJECTS += ../balloons-VALBAL/hootl/lasagna/PastaSim.cpp
HEADERS += ../balloons-VALBAL/hootl/lasagna/PastaSim.h



%.o: %.cpp $(HEADERS)
	$(CC) $(PATHS) $(CXXFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(PATHS) $(OBJECTS) -Wall $(LIBS) -o $@

diff: diff.o
	$(CC) diff.o -o difftest

clean:
	-rm -f *.o
	-rm -f $(TARGET)

open:
	vi -p traj.cpp config.h utils.h data.h data.cpp sim.h sim.cpp

opendiff:
	vi -p diff.h diff.cpp

install-adept:
	cd lib/adept-2.0.5 ; ./configure CXXFLAGS="-g1 -O3 -march=native" ; make -j2 ; sudo make install ; sudo echo "i totally couuld have just rm -rf the filesystem"

flight-code:
	cd .. ; ./utils/flight_code.sh

setup:
	cd .. ; sh utils/setup.sh
