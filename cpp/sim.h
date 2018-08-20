#ifndef SIM_H
#define SIM_H

#include <stdint.h>

#include "data.h"

template <class Float>
class PressureSource {
public:
	virtual Float get_pressure(int) = 0;
};

template <class Float>
class WindSource {
public:
	virtual wind_vector<Float> get_wind(int t, float lat, float lon, float pres) = 0;
};

template <class Float>
class Simulation {
public:
	Simulation(PressureSource<Float>& s, WindSource<Float>& w, int i);
	PressureSource<Float>& pressure;
	WindSource<Float>& winds;

	void run(int, float, float);

	int cur_file;

	const int dt = 60*10;
	const int tmax = 60*60*100;

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
	float *alts;
};

template <class Float>
class NearestNeighborWind : public WindSource<Float> {
public:
	wind_vector<Float> get_wind(int t, Float lat, Float lon, Float pres);
	int cur_file = 0;
};

#endif