#ifndef SIM_H
#define SIM_H

#include <random>
#include <stdint.h>

#include <adept.h>
using adept::adouble;

#include "data.h"	
#include "utils.h"
#include "opt.h"

#include "../ignored/balloons-VALBAL/src/LasagnaController.h"
#include "../ignored/balloons-VALBAL/hootl/lasagna/PastaSim.h"

/**
 * Basic pressure to altitude conversion and back
 * Meters, Pascals
 */
float p2alt(float p);
float alt2p(float alt);


template <class Float>
class Integrator {
public: 
	virtual void integrate(sim_state<Float>&, wind_vector<Float>&) = 0;
};

template <class Float>
class EulerInt : public Integrator<Float> {
public: 
	EulerInt(int dt) : dt(dt) {idlat = dt / (2 * M_PI * 6371008 / 360.);};
	EulerInt(){idlat = dt / (2 * M_PI * 6371008 / 360.);};
	void integrate(sim_state<Float>&, wind_vector<Float>&);
	int dt = 60*10;
	float idlat;
};

template <class Float>
class PressureSource {
public:
	virtual Float get_pressure(int, float, float) = 0;
};

template <class Float>
class WindSource {
public:
	virtual wind_vector<Float> get_wind(int, Float, Float, Float) = 0;
};


template <class Float>
class LinInterpWind : public WindSource<Float> {
public:
	LinInterpWind() : random_gen((std::random_device())()), normal(0,1) {};
	wind_vector<Float> get_wind(int, Float, Float, Float);
	int cur_file = 0;
	float sigma = 0;
    std::mt19937 random_gen;
    std::normal_distribution<> normal;
};


template <class Float>
class Simulation {
public:
	LinInterpWind<Float> wind_default;
	NoOp<Float> obj_default;
	EulerInt<Float> intg_default;
	Simulation(PressureSource<Float>& s, WindSource<Float>& w, ObjectiveFn<Float>& o, int i=-1);
	Simulation(PressureSource<Float>& s, int i=-1);
	PressureSource<Float>& pressure;
	WindSource<Float>& wind;
	Integrator<Float>& intg;
	ObjectiveFn<Float>& objfn;
	bool calc_obj = false;

	vec2<Float> run(int, Float, Float);

	int cur_file = 0;
	//const int tmax = 60*60*103;
	int tmax = 60*60*50;

	bool save_to_file = false;
	FILE *file;
private: 
	void init(int);
};

template <class Float>
class PressureTable : public PressureSource<Float> {
public:
	PressureTable(const char *);
	Float get_pressure(int, float, float);

	uint32_t t0;
	uint32_t dt;
	uint32_t n;
	float *alts;
};

template <class Float>
class WaypointController : public PressureSource<Float> {
public:
	WaypointController(int, int, Float *);
	Float get_pressure(int, float, float);

	int t0;
	int dt;
	Float *alts;
};


template <class Float>
class GreedySearch : public PressureSource<float> {
public:
	GreedySearch(WindSource<Float>& w, Integrator<Float>& i, ObjectiveFn<Float>& o, vec2<float>r) : wind(w), intg(i), objfn(o), range{alt2p(r.a),alt2p(r.b)} {};
	Float get_pressure(int, Float, Float);
	WindSource<Float>& wind;
	Integrator<Float>& intg;
	ObjectiveFn<Float>& objfn;
	vec2<Float> range;
	bool is_first_run = true;

	int N_levels = 100;

};

/**
 * Simulator for a valbal using a LasagnaController and a PastaSim from the balloons-VALBAL repo.
 */
template <class Float>
class LasSim : public PressureSource<Float> {
public: 
	LasSim(int seed);
	LasSim(int seed, float h);
	LasSim(int seed, float h, float l);
	Float get_pressure(int, float, float);
	const float freq = 1/60.;
	LasagnaController las;
	PastaSim sim;
	int t_last;
	bool is_first_run = true;
};

#endif
