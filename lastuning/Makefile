CC=g++
CFLAGS= -std=c++17 -D JOHNSIM -g -Ofast -I .
_PATHS = src lib/spa hootl/lasagna
_SRCF = src/Utils src/LasagnaController hootl/lasagna/PastaSim lib/spa/spa 
SRCF = $(addprefix ../ignored/balloons-VALBAL/, $(addsuffix .cpp, $(_SRCF)))
PATHS = $(addprefix -I ../ignored/balloons-VALBAL/, $(_PATHS))

all: build

build:
	$(CC) $(CFLAGS) $(PATHS) -o lasTune lasTuning.cpp $(SRCF)

run:
	./lasTune