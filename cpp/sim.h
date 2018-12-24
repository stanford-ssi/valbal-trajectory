#ifndef SIM_H
#define SIM_H

#include <random>
#include <stdint.h>

#include <adept.h>
using adept::adouble;

#include "data.h"	
#include "utils.h"
#include "objectives.h"
#include "trajtypes.h"
#include "opt.h"

#include "../ignored/balloons-VALBAL/src/LasagnaController.h"
#include "../ignored/balloons-VALBAL/hootl/lasagna/PastaSim.h"

/**
 * Basic pressure to altitude conversion and back
 * Meters, Pascals
 */
template<class Float>
Float p2alt(Float p){
	return (1.0-(pow((p/101350.0),0.190284)))*145366.45*0.3048;
}

template<class Float>
Float alt2p(Float alt){
	return pow(-((alt/145366.45/0.3048)-1.0),1.0/0.190284)*101350.0;
}


/**
 * Integrator object base class, which takes in the state and the wind, and propagates 
 * the simulation forward in time.
 */
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
	int dt = 10*60;
	float idlat;
};

template <class Float>
class EulerIntBal : public Integrator<Float> {
public: 
	EulerIntBal(int dt) : dt(dt) {idlat = dt / (2 * M_PI * 6371008 / 360.);};
	EulerIntBal(){idlat = dt / (2 * M_PI * 6371008 / 360.);};
	void integrate(sim_state<Float>&, wind_vector<Float>&);
	int dt = 60*10;
	float idlat;
};

/**
 * Pressure source base class. Takes in the state of the simulation and returns a pressure 
 * (which corrisponds to an altitude). The simulation runs in pressure and not altitude because
 * wind data is given in terms of pressure.
 */
template <class Float>
class PressureSource {
public:
	virtual void get_pressure(sim_state<Float>&) = 0;
};

/**
 * Wind source base class. Takes in the state of the simulation loads wind data, interpolates it, and
 * returns the wind vector for the given state
 */
template <class Float>
class WindSource {
public:
	virtual wind_vector<Float> get_wind(sim_state<Float>&) = 0;
};

/**
 * Basic wind source, uses simple linear interpolation of the data to return a wind
 */
template <class Float>
class LinInterpWind : public WindSource<Float> {
public:
	LinInterpWind(DataHandler& data) : dat(data), random_gen((std::random_device())()), normal(0,1) {};
	wind_vector<Float> get_wind(sim_state<Float>&);
	DataHandler& dat;
	int cur_file = 0;
	float sigma = 0;
    std::mt19937 random_gen;
    std::normal_distribution<> normal;
};

/**
 * Simulation base 
 */
template <class Float>
class Simulation {
public:
	LinInterpWind<Float> wind_default;
	NoOp<Float> obj_default;
	EulerInt<Float> intg_default;
	DataHandler dat_default;
	Simulation(PressureSource<Float>& s, WindSource<Float>& w, ObjectiveFn<Float>& o, Integrator<Float>& in, int i=-1);
	Simulation(PressureSource<Float>& s, WindSource<Float>& w, ObjectiveFn<Float>& o, int i=-1);
	Simulation(PressureSource<Float>& s, DataHandler& d, int i=-1);
	PressureSource<Float>& pressure;
	WindSource<Float>& wind;
	Integrator<Float>& intg;
	ObjectiveFn<Float>& objfn;
	bool calc_obj = false;

	sim_state<Float> run(int, Float, Float);
	void run(sim_state<Float>& state);
	
	int cur_file = 0;
	//const int tmax = 60*60*103;
	int tmax = 60*60*100;

	bool save_to_file = false;
	FILE *file;
private: 
	void init(int);
};

template <class Float>
class PressureTable : public PressureSource<Float> {
public:
	PressureTable(const char *);
	void get_pressure(sim_state<Float>&);

	uint32_t t0;
	uint32_t dt;
	uint32_t n;
	float *alts;
};

template <class Float>
class WaypointController : public PressureSource<Float> {
public:
	WaypointController(int, int, Float *);
	void get_pressure(sim_state<Float>&);

	int t0;
	int dt;
	Float *alts;
};


template <class Float>
class GreedySearch : public PressureSource<float> {
public:
	GreedySearch(WindSource<Float>& w, Integrator<Float>& i, ObjectiveFn<Float>& o, vec2<float>r) : wind(w), intg(i), objfn(o), range{alt2p(r.a),alt2p(r.b)} {};
	void get_pressure(sim_state<Float>&);
	WindSource<Float>& wind;
	Integrator<Float>& intg;
	ObjectiveFn<Float>& objfn;
	vec2<Float> range;
	bool is_first_run = true;

	int N_levels = 100;
};


template<class Float>
class ParameterServer {
public:
	virtual ctrl_cmd<Float> get_param(sim_state<Float>&) = 0;
	virtual double apply_gradients(StepRule&) = 0;
};

/* I hate C++. --Joan */
template <typename> struct tag {};

template<class Float>
class TemporalParameters : public ParameterServer<Float> {
public:
	TemporalParameters(int t0_, int dt_, int T_, double d_h, double d_t);
	~TemporalParameters();
	ctrl_cmd<Float> get_param(sim_state<Float>&);
	double apply_gradients(StepRule&);
	template <class FFloat> TemporalParameters& operator = (const TemporalParameters<FFloat>& tp){
		for(int i=0; i < int(T/dt); i++){
			this->cmds[i].h = VAL(tp.cmds[i].h);
			this->cmds[i].tol = VAL(tp.cmds[i].tol);
		}
		return *this;
	}
	ctrl_cmd<Float> *cmds;
	int t0;
	int dt;
	int T;
	int N;
private:
	double default_h;
	double default_tol;
	double apply_gradients(StepRule&, tag<TemporalParameters<float>>);
	double apply_gradients(StepRule&, tag<TemporalParameters<adouble>>);
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
	LasSim(int seed, float h, float l, TemporalParameters<float>& cmds);
	void get_pressure(sim_state<Float>& state);
	void get_altitude(sim_state<Float>& state);
	const float freq = 1/60.;
	LasagnaController las;
	PastaSim sim;
	int t_last = 0;
	int dt = 0;
	bool changing_cmds = false;
	LasagnaController::Constants lasconst;
	TemporalParameters<float> cmds_defualt;
	TemporalParameters<float>& cmds;
private: 
	void evolve(sim_state<Float>& state); 
};

template <class Float>
class StochasticControllerApprox : public PressureSource<Float> {
public: 
	StochasticControllerApprox(ParameterServer<Float>& ps_, int seed);
	void get_pressure(sim_state<Float>&);
	LasSim<float> las_sim;
	float h_mid;
	float tol0;
	ParameterServer<Float>& params;
};

template<class Float>
class cmd_tree {
public:
	cmd_tree() : upper(0), lower(0) {
		printf("Initializing cmd tree!\n");
		for (int i=0; i<D; i++) {
			mins[i] = 1e20;
			maxs[i] = -1e20;
		}
	};
	~cmd_tree() {
		delete upper;
		delete lower;
	}
	template <class FFloat> cmd_tree& operator = (const cmd_tree<FFloat>& tp){
		if (tp.upper != 0) {
			cmd_tree<Float> *up = new cmd_tree<Float>;
			*up = *tp.upper;
			this->upper = up;
		} else {
			this->upper = 0;
		}
		if (tp.lower != 0) {
			cmd_tree<Float> *lo = new cmd_tree<Float>;
			*lo = *tp.lower;
			this->lower = lo;
		} else {
			this->lower = 0;
		}
		this->a = tp.a;
		this->b = tp.b;
		this->c = tp.c;
		this->cmd.h = VAL(tp.cmd.h);
		this->cmd.tol = VAL(tp.cmd.tol);
		return *this;
	}

	void save_to_file(const char *fname);
	void _save_to_file(FILE *f);

	ctrl_cmd<Float>& get_cmd(sim_state<Float>&);
	void gradients_and_split(StepRule&);

	cmd_tree<Float> *upper;
	cmd_tree<Float> *lower;

	ctrl_cmd<Float> cmd;

	static const int D = 2;

	vector<array<float, D>> requests;
	array<float, D> mins;
	array<float, D> maxs;

	// a*lat + b*lon >= c --> upper
	double a;
	double b;
	double c;
};

template<class Float>
class SpatiotemporalParameters : public ParameterServer<Float> {
public:
	SpatiotemporalParameters(int t0_, int dt_, int T_, double d_h, double d_t);
	template <class FFloat> SpatiotemporalParameters& operator = (const SpatiotemporalParameters<FFloat>& tp){
		for(int i=0; i < int(T/dt); i++){
			this->cmds[i] = tp.cmds[i];
		}
		return *this;
	}
	~SpatiotemporalParameters();
	ctrl_cmd<Float> get_param(sim_state<Float>&);
	double apply_gradients(StepRule&);

	int t0;
	int dt;
	int T;
	cmd_tree<Float> *cmds;
	double default_h;
	double default_tol;
	double apply_gradients(StepRule&, tag<SpatiotemporalParameters<float>>);
	double apply_gradients(StepRule&, tag<SpatiotemporalParameters<adouble>>);
};

#endif
