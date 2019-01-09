#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "LasagnaController.h"
#include "header.h"
#include "PastaSim.h"

#define CONTROLLER LasagnaController
#define LAS LasagnaController

//#define WRITE

using namespace std;


CONTROLLER::Constants getDefaultConstants();
CONTROLLER::Constants DefaultConstants();

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}


int main(int argc, char *argv[])
{


	LasagnaController::Constants c;
	char* p;
	p = getCmdOption(argv, argv + argc, "--tol");
	if(p) c.tolerance = atof(p);
	p = getCmdOption(argv, argv + argc, "--gain");
	if(p) c.gain = atof(p);
	p = getCmdOption(argv, argv + argc, "--damp");
	if(p) c.damping = atof(p);
	
	int days = 100;
	p = getCmdOption(argv, argv + argc, "--days");
	if(p) days = atoi(p);

	p = getCmdOption(argv, argv + argc, "--seed");
	int seed = 0;
	if(p) seed = atoi(p);
	PastaSim sim = p ? PastaSim(seed) : PastaSim();

	bool write = cmdOptionExists(argv, argv+argc, "-w");

	//printf("%f,%d\n",tol,argc);
	const float freq = 1;
	fstream o ("output.bin", std::fstream::out | std::fstream::binary);
	
	sim.conf.nightfall = 0;
	sim.conf.freq = freq;
	sim.h = 13500;
	sim.l = 0;
	sim.conf.k_drag = 0.006;
	CONTROLLER las(freq);
	las.updateConstants(c);
	miniframe data;
	float v_cmd = 0;
	int dur = 60*60*24*days*freq;
	int act_sum = 0;
	float bal_sum = 0;
	for(int i = 0; i < 60*60*20*4; i++){
		CONTROLLER::Input input;
		input.h_rel = sim.h;
		input.h_abs = sim.h;		
		las.update(input);
	}
	for(int i = 0; i < dur; i++){
		CONTROLLER::Input input;
		float h = sim.evolve(double(las.getAction()));
		float bal = (las.getAction() > 0)*las.getAction()/1000.*sim.conf.bal_dldt;
		bal_sum += bal;
		input.h_rel = h;
		input.h_abs = h;
		input.dldt_ext = 4*sim.sunset_dldt;
		las.update(input);
		CONTROLLER::State state = las.getState();
		act_sum += state.action;
		if(write){
			if(i%(int(ceil(freq))) == 0){
				float buf[1] = {sim.h};
				o.write((char*)&buf, sizeof(buf));
				//o.write((char*)&act_sum,sizeof(act_sum));
			}
		}
		//if(i%int(freq*60*60*24) == 0) printf("%f ", i/freq/60/60);
	}
	printf("%f\n",bal_sum/dur*60*60);
}
