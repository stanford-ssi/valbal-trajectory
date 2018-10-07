#include <stdio.h>
#include <assert.h>

#include <adept.h>
using adept::adouble;

#include "data.h"
#include "sim.h"
#include "utils.h"

#define STORE_ALTITUDE

using namespace std;

void simpleSim(){
	//PressureTable<float> pres(get_data_path("flights/ssi71_inp.bin"));
	TIMEIT("Running simple sims",
		const int N = 10000;
    	Scheduler<float> sched(-2, N);
		for (int i=0; i<N; i++) {
			sched.add([i]() {
				LasSim<float> pres(std::time(0)+i,1000,0.1);
				Simulation<float> sim(pres, i);
				sim.wind_default.sigma = 2;
				sim.tmax = 60*60*100;
				vec2<float> f = sim.run(1536994800, 36.849014, -121.432913+360);
				return f.a;
			});
		}
		sched.finish();
	)
	for (int i=0; i<N; i++) {
		assert(sched.results[i] != 0);
	}
}

/* Lol ignore these functions it's joan being lazy */
/* lol gdi joan --john */
void simpleSim63(){
	PressureTable<float> pres(get_data_path("flights/ssi63_inp.bin"));
	TIMEIT("Running simple sims",
		const int N = 200;
		Scheduler<float> sched(-1, N);
		for (int i=0; i<N; i++) {
			sched.add([&pres, i]() {
				Simulation<float> sim(pres, i);
				sim.wind_default.sigma = 1;
				vec2<float> f = sim.run(1512889140+7*60*60, 37.251022, -122.03919+360);
				return f.a;
			});
		}
		sched.finish();
	)
	for (int i=0; i<N; i++) {
		assert(sched.results[i] != 0);
	}
}

void simpleSim67(){
	PressureTable<float> pres(get_data_path("flights/ssi67_inp.bin"));
	TIMEIT("Running simple sims",
		const int N = 200;
    	Scheduler<float> sched(-1, N);
		for (int i=0; i<N; i++) {
			sched.add([&pres, i]() {
				Simulation<float> sim(pres, i+1);
				if (i != 0) sim.wind_default.sigma = 1;
				vec2<float> f = sim.run(1526169840, 35.714558, -119.94677+360);
				return f.a;
			});
		}
		sched.finish();
	)
	for (int i=0; i<N; i++) {
		assert(sched.results[i] != 0);
	}
}


void ssi71Sims(){
	// for playing around with ssi71 data
	PressureTable<float> pres1(get_data_path("flights/ssi71_inp.bin"));
	Simulation<float> sim(pres1, 1);
	vec2<float> f = sim.run(1529007840, 42.46636, -68.021255+360);
	printf("SIMDONE %f %f\n", f.a, f.b);
	for(int i = 2; i < 50; i++){
		LasSim<float> pres2(13500.);
		pres2.sim.l = 0.0;
		pres2.sim.conf.gtime.year = 2018;
		pres2.sim.conf.gtime.month = 6;
		pres2.sim.conf.gtime.day = 14;
		pres2.sim.conf.gtime.hour = 20;
		pres2.sim.conf.gtime.minute = 34;
		pres2.sim.conf.gtime.second = 0;
		//LasagnaController::Constants con; con.h_cmd = 
		//pres2.las.update
		Simulation<float> sim(pres2, i);
		vec2<float> f = sim.run(1529007840, 42.46636, -68.021255+360);
		printf("SIMDONE %f %f\n", f.a, f.b);
	}

}

void gradientsStuff(){
	int t0 = 1536562800;
	int dt = 3600*5;
	const int N_W = 21;
	const float LR = 10; (void)LR;
	double waypoints_val[N_W]; for (int i=0; i<N_W; i++) waypoints_val[i] = alt2p(18000);

	adept::Stack stack;
	for (int it=0; it<1000; it++) {
		clock_t timer0 = clock();
		adouble waypoints[N_W];
		for (int i=0; i<N_W; i++) {
			waypoints[i] = min(double(alt2p(5000.)), max(double(alt2p(19000.)), waypoints_val[i]));
		}
		stack.new_recording();
	
		WaypointController<adouble> pres(t0, dt, waypoints);
		LinInterpWind<adouble> wind;
		//FinalLongitude<adouble> obj;

		MinDistanceToPoint<adouble> obj(31.753952, -77.081033+360);
		Simulation<adouble> sim(pres, wind, obj,it+1);
		sim.tmax=60*60*100;
		vec2<adouble> end = sim.run(pres.t0, 40.785741, -73.410338+360);
		adouble cost = -obj.getObjective()*111195;

		cost.set_gradient(1.0);
		stack.compute_adjoint();
		for (int i=0; i<N_W; i++) {
			waypoints_val[i] += LR * waypoints[i].get_gradient();
			printf("%.1f, ",p2alt(waypoints_val[i])/1000);
		} 
		printf("\n");
		float dt = (clock() - timer0)/((double)CLOCKS_PER_SEC)*1000;
		printf("Took %.2f ms, got %f\n", dt, VAL(cost)/1e6);
	}
}

void searchStuff(){
	clock_t timer0 = clock();
	//MinDistanceToPoint<float> obj(59.916193, 30.325234+360);
	FinalLongitude<float> obj;
	LinInterpWind<float> wind;
	EulerInt<float> intg;
	GreedySearch<float> pres(wind,intg,obj,vec2<float>{10000,18000});
	Simulation<float> sim(pres, wind, obj,0);
	sim.tmax=60*60*100;
	sim.run(1536994800, 36.95854187, -121.84505463+360);
	float dt = (clock() - timer0)/((double)CLOCKS_PER_SEC)*1000;
	printf("Took %.2f ms, got %f\n", dt, obj.getObjective()/1e6);
}


void stocasticGradients(){
	int t0 = 1536562800;
	int dt = 3600*5;
	const int N_W = 21;
	const float LR = 10; (void)LR;
	ctrl_cmd<float> cmd;
	cmd.h = 13000;
	cmd.tol = 750; 
	ctrl_cmd<float> cmds[N_W]; 
	for (int i=0; i<N_W/4; i++) cmds[i] = cmd;
	cmd.h = 14000;
	cmd.tol = 350;
	for (int i=N_W/4; i<N_W/2; i++) cmds[i] = cmd;
	cmd.h = 12000;
	cmd.tol = 1000;
	for (int i=N_W/2; i<3*N_W/4; i++) cmds[i] = cmd;
	cmd.h = 14500;
	cmd.tol = 500;
	cmds[3*N_W/4] = cmd;
	cmd.h = 17000;
	cmd.tol = 500;
	for (int i=3*N_W/4+1; i<N_W; i++) cmds[i] = cmd;


	for (int it=0; it<500; it++) {
		clock_t timer0 = clock();
		StocasticControllerApprox<float> controller(t0, dt, cmds, std::time(0)+it*1000);
		LinInterpWind<float> wind;
		NoOp<float> obj;
		Simulation<float> sim(controller, wind, obj, it);
		sim.tmax=60*60*100;
		sim.run(controller.t0, 40.785741, -73.410338+360);
		float dt = (clock() - timer0)/((double)CLOCKS_PER_SEC)*1000;
		printf("Took %.2f ms\n",dt);
	}
}

int main() {
	printf("ValBal trajectory optimization.\n");
	printf("This program is highly cunning and, on the face of it, not entirely wrong.\n");

	load_data(get_data_path("proc/gfs_anl_0deg5"), 1500000000,1600000000);
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
	stocasticGradients();
}
