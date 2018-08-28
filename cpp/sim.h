#ifndef SIM_H
#define SIM_H

#include <stdint.h>

#include <adept.h>
using adept::adouble;

#include "data.h"

#include "../ignored/balloons-VALBAL/src/LasagnaController.h"
#include "../ignored/balloons-VALBAL/hootl/lasagna/PastaSim.h"

inline float VAL(float f) { return f; }
inline float VAL(adouble f) { return f.value(); }

template <class Float>
class PressureSource {
public:
	virtual Float get_pressure(int) = 0;
};

template <class Float>
class WindSource {
public:
	virtual wind_vector<Float> get_wind(int t, Float lat, Float lon, Float pres) = 0;
};

template <class Float>
class Simulation {
public:
	Simulation(PressureSource<Float>& s, WindSource<Float>& w, int i);
	PressureSource<Float>& pressure;
	WindSource<Float>& winds;

	vec2<Float> run(int, Float, Float);

	int cur_file;

	const int dt = 60*10;
	const int tmax = 60*60*50;

	FILE *file;
};

template <class Float>
class PressureTable : public PressureSource<Float> {
public:
	PressureTable(const char *);
	Float get_pressure(int);

	uint32_t t0;
	uint32_t dt;
	uint32_t n;
	float *alts;
};

template <class Float>
class WaypointController : public PressureSource<Float> {
public:
	WaypointController(int, int, Float *);
	Float get_pressure(int);

	int t0;
	int dt;
	Float *alts;
};

/**
 * Simulator for a valbal using a LasagnaController and a PastaSim from the balloons-VALBAL repo.
 */
template <class Float>
class LasSim : public PressureSource<Float> {
public: 
	LasSim(float start_h) {sim.h = start_h;}; 
	Float get_pressure(int);
	LasagnaController las;
	PastaSim sim;
	int t_last;
	bool is_first_run = true;
};

template <class Float>
class NearestNeighborWind : public WindSource<Float> {
public:
	wind_vector<Float> get_wind(int t, Float lat, Float lon, Float pres);
	int cur_file = 0;
};


/**
 * Basic pressure to altitude conversion and back
 * Meters, Pascals
 */
float p2alt(float p);
float alt2p(float alt);


#endif
