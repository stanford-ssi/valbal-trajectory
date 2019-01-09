#include <stdio.h>
#include <assert.h>

#include <adept.h>
using adept::adouble;

#include "data.h"
#include "sim.h"
#include "utils.h"
#include "planners.h"
#include "trajtypes.h"
#include <random>
#include <dirent.h>
#define STORE_ALTITUDE

using namespace std;

void demo() {
	DataHandler data;
	data.load_data(get_data_path("/proc/gfs_anl_0deg5"), 1500000000,1600000000);
	int t0 = 1540755060;
	float lat0 = 36.84;
	float lon0 =  -121.43 + 360;

	/* Fly at constant 14 km. This is the most simple controller. Other controllers
     * include simulating the Lasagna controller, or having a differentiable controller
     * that can be optimized. */
	float pressures[150];
	for (int i=0; i<150; i++) pressures[i] = alt2p(14000);
	WaypointController<float> controller(t0, 3600, pressures);

    TIMEIT("Running simple sims",
        const int N = 5000000;
		/* This is just to run it faster, on multiplee threads. The scheduler isn't
         * strictly necessary. */
        Scheduler<float> sched(-2, N);
        for (int i=0; i<N; i++) {
            sched.add([&controller, i, t0, lat0, lon0, &data]() {
                Simulation<float> sim(controller, data, -1);
                sim.wind_default.sigma = 0;
				sim.tmax = 60*60*60;
                return sim.run(t0, lat0, lon0).lon;
            });
        }
        sched.finish();
    )
    for (int i=0; i<N; i++) {
        assert(sched.results[i] != 0);
    }

}

void stochasticGradients(){
	//const char* db = get_data_path("/proc/gfs_pred_0deg5/20181129_12/");
	const char* db = get_data_path("/proc/gfs_anl_0deg5/");
	sim_state<float> state0;
	state0.lat = 36.84;
	state0.lon = -121.43 + 360;
	state0.t = 1543492801;
	state0.bal = 4.5;
	StochasticMPC controller(db,state0);
	controller.conf.opt_sign = 1;
	controller.run();
}

void spatialGradients(){
	//const char* db = get_data_path("/proc/gfs_pred_0deg5/20181129_12/");
	const char* db = get_data_path("/proc/gfs_anl_0deg5/");
	sim_state<float> state0;
	state0.lat = 36.84;
	state0.lon = -121.43 + 360;
	state0.t = 1543492801;
	state0.t = 1541041200;
	state0.bal = 4.5;
	SpatialPlanner controller(db,state0);
	controller.conf.opt_sign = 1;
	controller.run();
}

void evaluator(){
	DataHandler anldata;
	anldata.load_data(get_data_path("proc/gfs_anl_0deg5/"),1500000000,1600000000);
	sim_state<float> state;
	state.lat = 36.84;
	state.lon = -121.43 + 360;
	int t0= 1543492801;
	state.t = t0;
	state.bal = 4.5;
	state.p = alt2p(13000);
	float lift = 0; 
	const int command_intval = 60*60;

	sim_state<float> statecpy = state;
	TemporalParameters<float> cmds(statecpy.t, 60*60*6, 60*60*150, 13000, 2000);
	LasSim<float> alt_sim(0,p2alt(statecpy.p),lift,cmds);
	EulerIntBal<float> in;
	alt_sim.dt = in.dt; 
	NoOp<float> obj;
	LinInterpWind<float> anlwind(anldata);
	Simulation<float> sim(alt_sim,anlwind,obj,in,0);
	sim.tmax = 60*60*150;
	sim.run(statecpy);


	int i = 0;
	while(state.bal > 0){
		char recentdir[PATH_MAX];
		getRecentDir(recentdir,get_data_path("proc/gfs_pred_0deg5/"),state.t);
		printf("%s\n",recentdir);
		StochasticMPC controller(recentdir,state);
		//controller.conf.n_iters = 20;
		controller.conf.n_samples = 50;
		controller.conf.write_files = true;
		controller.conf.opt_sign = -1;
		controller.conf.fname_offset = 100+i*2*controller.conf.n_samples;
		TemporalParameters<float> cmds = controller.run();
		LasSim<float> alt_sim(i+1,p2alt(state.p),lift,cmds);
		EulerIntBal<float> in;
		alt_sim.dt = in.dt; 
		NoOp<float> obj;
		LinInterpWind<float> anlwind(anldata);
		Simulation<float> sim(alt_sim,anlwind,obj,in,i+1);
		sim.tmax = command_intval;
		sim.run(state);
		lift = alt_sim.sim.l;
		printf("###############################################################\n[real state] %d lat:%f, lon:%f, bal:%f, alt:%f\n",i,state.lat,state.lon,state.bal,p2alt(state.p));
		i++;
	}

}

sim_state<adept::adouble> mkstate(float lat, float lon) {
	sim_state<adept::adouble> out;
	out.lat = lat;
	out.lon = lon;
	return out;
}

void print(ctrl_cmd<adept::adouble> s) {
	printf("got command (%f, %f)\n", VAL(s.h), VAL(s.tol));
}

void test_spatial() {

	adept::Stack stack;
	cmd_tree<adept::adouble> tree;
	tree.cmd.h = 14000;
	tree.cmd.tol = 1500;

	auto s = mkstate(42, 2);
	print(tree.get_cmd(s));
	s = mkstate(44, 2);
	print(tree.get_cmd(s));
	s = mkstate(44.1, 3);
	print(tree.get_cmd(s));

	adept::adouble obj = 3.1;
	stack.new_recording();	
	obj.set_gradient(1.0);

	stack.compute_adjoint();
	GradStep step;
	tree.gradients_and_split(step);

}

int main() {
	printf("ValBal trajectory optimization.\n");
	printf("This program is highly cunning and, on the face of it, not entirely wrong.\n");


	//load_data(get_data_path("/proc/gfs_pred_0deg5/20181021_12"), 1500000000,1600000000);
	//load_data("../ignored/proc/euro_anl", 1500000000,1600000000);
	//load_data("../ignored/proc/euro_fc", 1500000000,1600000000);
	//load_data("../proc", 1500000000,1600000000);

	/*wind_data *sample = get_data_at_point(files+0, {42, 120});
	printf("(u,v) %hd %hd\n", sample->data[4][0], sample->data[4][1]);

	point base = get_base_neighbor(69.5, 60.6);
	point near = get_nearest_neighbor(69.5, 60.9);
	printf("(%d,%d) (%d, %d)\n", base.lat, base.lon, near.lat, near.lon);
	*/
	//simpleSim67();
	//ssi71Sims();
	//gradientsStuff();
	//MLestimation();
	stochasticGradients();
	//test_spatial();
	//spatialGradients();

	//evaluator();
	//saveSpaceshot();
	//stocasticGradients();
}
