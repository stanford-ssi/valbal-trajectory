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
	const char* db = get_data_path("/proc/gfs_pred_0deg5/20181129_12/");
	sim_state<float> state0;
	state0.lat = 36.84;
	state0.lon = -121.43 + 360;
	state0.t = 1543492800;
	StocasticMPC controller(db,state0);
	//controller.conf.opt_sign = -1;
	controller.run();
}

void evaluator(){
	DataHandler anldata;
	anldata.load_data(get_data_path("proc/gfs_anl_0deg5/"),1500000000,1600000000);
	sim_state<float> state;
	state.lat = 36.84;
	state.lon = -121.43 + 360;
	state.t = 1543492801;
	state.bal = 4.5;



	char recentdir[PATH_MAX];
	getRecentDir(recentdir,get_data_path("proc/gfs_pred_0deg5/"),state.t);
	printf("%s\n",recentdir);
	StocasticMPC controller(recentdir,state);
	controller.run();
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
	//stochasticGradients();
	evaluator();
	//saveSpaceshot();
	//stocasticGradients();
}
