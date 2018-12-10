#include <stdio.h>
#include <assert.h>

#include <adept.h>
using adept::adouble;

#include "data.h"
#include "sim.h"
#include "utils.h"
#include <random>
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
	DataHandler data;
	data.load_data(get_data_path("../ignored/proc/gfs_pred_0deg5/20181129_12/"), 1500000000,1600000000);
    //load_data(get_data_path("../ignored/proc/gfs_anl_0deg5"), 1500000000,1600000000);
    //int t0 = 1537196400;
    int t0 = 1543492800;
	float lat0 = 36.84;
	float lon0 =  -121.43 + 360;
	int dt = 3600*6;
	const int N_RUNS = 50;
	float LR_H = 200000;
	float LR_T = 20000;
	const float alpha = 0.995;
	const int N_IT = 200;
	adept::Stack stack;
	Optimizer opt(LR_H, LR_T);
	TemporalParameters<adouble> params(t0, dt, 120*dt, 14000, 2000);
	for (int it=0; it<N_IT; it++){
		LR_H = LR_H*alpha;
		clock_t timer0 = clock();
		adouble obj_sum = 0;
		stack.new_recording();
		adouble objectives[N_RUNS];
		float meanbal = 0;
		float meantime = 0;

		for (int run=0; run<N_RUNS; run++) {
			StochasticControllerApprox<adouble> controller(params, rand());
			LinInterpWind<adouble> wind(data);
            wind.sigma = 0;
			FinalLongitude<adouble> obj;
			EulerIntBal<adouble> in;
			int fname = -1;
			if (it == 0 || it == N_IT-1) { printf("saving!\n"); fname = run + N_RUNS*(it == 0); }
			Simulation<adouble> sim(controller, wind, obj, in, fname);
			sim.tmax=60*60*120;
			sim_state<adouble> sf = sim.run(t0, lat0, lon0);
			meanbal += VAL(sf.bal);
			meantime += (sf.t - t0);
			objectives[run] = obj.getObjective();
		}
		meanbal /= N_RUNS;
		meantime /= N_RUNS;

		for (int run=0; run<N_RUNS; run++) obj_sum += objectives[run];
		obj_sum = obj_sum/((float)N_RUNS);
		obj_sum.set_gradient(1.0);
		stack.compute_adjoint();
		params.apply_gradients(opt);

		float dt = (clock() - timer0)/((double)CLOCKS_PER_SEC)*1000;
		printf("Took %.2f ms, obj: %f, bal: %f, days: %f\n",dt,VAL(obj_sum), meanbal, meantime/86400.);
	}
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
	//saveSpaceshot();
	//stocasticGradients();
}
