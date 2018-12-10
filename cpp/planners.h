#ifndef PLANNERS_H
#define PLANNERS_H


#include "sim.h"
#include "objectives.h"

class StocasticMPC {
public: 
	typedef struct{
		double lat0 = 36.84;
		double lon0 = -121.43 + 360;
		double tmax = 120;
		double simdt = 600*6;
		int n_samples = 50;
		int n_iters = 300;
	} Config;
	StocasticMPC(const char* intput_db, int t0){};
	DataHandler data;
	adept::Stack stack;
	StepRule opt;
};


#endif