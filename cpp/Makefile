TARGET = traj
LIBS = -lm -ladept -lpthread -ljsoncpp
CC = g++ -flto
#-fprofile-generate
#-fprofile-use -fprofile-correction
CXXFLAGS = -g -Werror -Ofast -march=native -fno-strict-aliasing -Wall -ffast-math -std=c++14

PATHS = -I ../ignored/balloons-VALBAL/lib/spa
PATHS += -I ../ignored/balloons-VALBAL/src
PATHS += -I /usr/local/include


.PHONY: default all clean

default: $(TARGET)
all: default
run: default
	./traj

OBJECTS = $(patsubst %.cpp, %.o, $(filter-out extra.cpp, $(wildcard *.cpp)))
HEADERS = $(wildcard *.h)

#Jank additions to compile lasagna and Pastasim. Must run utils/flight_code.sh first
OBJECTS += ../ignored/balloons-VALBAL/src/LasagnaController.cpp
HEADERS += ../ignored/balloons-VALBAL/src/LasagnaController.h
OBJECTS += ../ignored/balloons-VALBAL/src/Utils.cpp
HEADERS += ../ignored/balloons-VALBAL/src/Utils.h
OBJECTS += ../ignored/balloons-VALBAL/lib/spa/spa.cpp
HEADERS += ../ignored/balloons-VALBAL/lib/spa/spa.h
OBJECTS += ../ignored/balloons-VALBAL/hootl/lasagna/PastaSim.cpp
HEADERS += ../ignored/balloons-VALBAL/hootl/lasagna/PastaSim.h



%.o: %.cpp $(HEADERS)
	$(CC) $(PATHS) $(CXXFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(PATHS) $(OBJECTS) -g -Wall $(LIBS) -o $@

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
