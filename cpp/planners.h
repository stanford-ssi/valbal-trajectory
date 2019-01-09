#ifndef PLANNERS_H
#define PLANNERS_H


#include "sim.h"
#include "objectives.h"


class StochasticMPC {
public: 
	typedef struct{
		sim_state<float> state0;
		int tmax = 120*60*60;				// max sim time in seconds
		int cmd_dt = 3600*6;				// time between comand waypoints
		int init_cmd_dt = 60*60*24;
		int n_samples = 50;					// Number of samples per batch
		int n_iters_max = 1000;				// Number of itterations max
		int n_iters_min = 10;				// Number of itterations min
		int n_starts = 50;					// Numer of random starts to try
		int opt_sign = 1;					// sign on the optimization
		int fname_offset = 0;				// offset on the file name numbers
		bool write_files = true;
	} Config;

	typedef struct {
		double lr_set  = 200000;
		double lr_tol  = 200000;
		double alpha   = 0.99999;
		double tol_min = 200;
		double tol_max = 3000;
		double set_min = 10000;
		double set_max = 16000;
		double obj_filter_corner = 0.03;
		double convergence_thresh = 0.1;
	} HyperParams;
	HyperParams hparams;
	FILE *file;
	bool save_to_file = 1;



	StochasticMPC(const char* input_db, sim_state<float> state0);
	TemporalParameters<float> run();
	DataHandler data;
	Config conf;
	//FinalLongitude<adouble> objfn;
};

class SpatialPlanner {
public: 
	typedef struct{
		sim_state<float> state0;
		int tmax = 120*60*60;				// max sim time in seconds
		int cmd_dt = 3600*6;				// time between comand waypoints
		int n_samples = 50;					// Number of samples per batch
		int n_iters = 300;					// Number of itterations
		int n_starts;						// Numer of random starts to try
		int opt_sign = -1;					// sign on the optimization
		int fname_offset = 0;				// offset on the file name numbers
		bool write_files = true;
	} Config;
	SpatialPlanner(const char* input_db, sim_state<float> state0);
	adept::Stack stack;
	SpatiotemporalParameters<float> run();
	DataHandler data;
	GradStep step;
	Config conf;
	//FinalLongitude<adouble> objfn;
};



#endif
