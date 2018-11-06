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
				sim.run(1536994800, 36.849014, -121.432913+360);
				return 0.;
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
				sim.run(1512889140+7*60*60, 37.251022, -122.03919+360);
				return 0.;
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
				sim.run(1526169840, 35.714558, -119.94677+360);
				return 0.;
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
	sim.run(1529007840, 42.46636, -68.021255+360);
	printf("SIMDONE");
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
		sim.run(1529007840, 42.46636, -68.021255+360);
		printf("SIMDONE\n");
	}

}

void gradientsStuff(){
	load_data(get_data_path("/proc/gfs_pred_0deg5/20181029_12"), 1500000000,1600000000);
	int t0 = 1540834728;
	int dt = 3600*5;
	const int N_W = 41;
	const float LR = 500000; (void)LR;
	const int N_ITER = 700;
  	std::default_random_engine gen(1);
  	std::uniform_real_distribution<float> init_alt(12500,19000);

	for(int j=0; j<100;j++){
		char path[100];
		snprintf(path, 100, "../ignored/sim/opt.%03d.bin", j);
		FILE* file = fopen(path, "wb");
		ensure(file != 0);
		ensure(setvbuf(file, 0, _IOFBF, 16384) == 0);
		double waypoints_val[N_W]; for (int i=0; i<N_W; i++) waypoints_val[i] = alt2p(init_alt(gen));
		adept::Stack stack;
		float cost_actual;
		for (int it=0; it<N_ITER; it++) {
			clock_t timer0 = clock();
			adouble waypoints[N_W];
			for (int i=0; i<N_W; i++) {
				waypoints[i] = min(double(alt2p(12500.)), max(double(alt2p(19000.)), waypoints_val[i]));
				//printf("%.1f, ",p2alt(VAL(waypoints[i]))/1000);
			}
			stack.new_recording();
		
			WaypointController<adouble> pres(t0, dt, waypoints);
			LinInterpWind<adouble> wind;
			FinalLongitude<adouble> obj;
			//MinDistanceToPoint<adouble> obj(31.753952, -77.081033+360);
			Simulation<adouble> sim(pres, wind, obj,N_ITER*j + it);
			sim.tmax=60*60*100;
			sim.run(pres.t0, 39.0152, -91.1433+360);
			adouble cost = obj.getObjective();

			cost.set_gradient(1.0);
			stack.compute_adjoint();
			for (int i=0; i<N_W; i++) {
				waypoints_val[i] += LR * waypoints[i].get_gradient();
			} 
			float dt = (clock() - timer0)/((double)CLOCKS_PER_SEC)*1000; (void)dt;
			cost_actual = VAL(cost);
			fwrite(&cost_actual, sizeof(float), 1, file);
		}
		printf("Final Cost: %f\n", cost_actual);
		fclose(file);
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
	load_data(get_data_path("/proc/gfs_pred_0deg5/20181101_00"), 1500000000,1600000000);
	int t0 = 1541045606;
	int dt = 3600*6;
	const int N_W = 26;
	const int N_RUNS = 50;
	const float LR_TOL = 10000;
	float LR_H = 2000000;
	const float alpha = 0.98;
	ctrl_cmd<double> cmd;
	cmd.h = 16000;
	cmd.tol = 1000; 
	ctrl_cmd<double> cmds_val[N_W]; 
	for (int i=0; i<N_W; i++) cmds_val[i] = cmd;
	//cmd.h = 13000;
	//for (int i=0; i<N_W/8; i++) cmds_val[i] = cmd;	
	adept::Stack stack;
	for (int it=0; it<300; it++){
		LR_H = LR_H*alpha;
		clock_t timer0 = clock();
		ctrl_cmd<adouble> cmds[N_W];
		for (int i=0; i<N_W; i++){
			cmds[i].h = cmds_val[i].h;
			cmds[i].tol = cmds_val[i].tol;	
		};
		adouble obj_sum = 0;
		stack.new_recording();
		int tf;
		for (int run=0; run<N_RUNS; run++) {
			StocasticControllerApprox<adouble> controller(t0, dt, cmds,run*1019);
			LinInterpWind<adouble> wind;
			FinalLongitude<adouble> obj;
			//MinDistanceToPoint<adouble> obj(40.378182, -3.958676+360);

			EulerIntBal<adouble> in;
			Simulation<adouble> sim(controller, wind, obj, in, run + N_RUNS*(it==0));
			sim.tmax=60*60*100;
			tf = sim.run(controller.t0,47.4289, -19.6931+360).t;
			obj_sum = obj.getObjective();
		}
		obj_sum = obj_sum/((float)N_RUNS);
		obj_sum.set_gradient(1.0);
		stack.compute_adjoint();

		for (int i=0; i<N_W; i++) {
			cmds_val[i].h += LR_H * cmds[i].h.get_gradient();
			cmds_val[i].h = min(16500., max(10000., cmds_val[i].h));
			//cmds_val[i].tol += LR_TOL * cmds[i].tol.get_gradient();
			//cmds_val[i].tol = min(2000., max(200., cmds_val[i].tol));
			printf("%.1f, ",cmds_val[i].h/1000);
		} printf("\n");
		for (int i=0; i<N_W; i++){printf("%.2f, ",LR_H * cmds[i].h.get_gradient());} printf("\n");
		for (int i=0; i<N_W; i++){printf("%.2f, ",cmds_val[i].tol/1000);} printf("\n");
		for (int i=0; i<N_W; i++){printf("%.2f, ",LR_TOL*cmds[i].tol.get_gradient());} printf("\n");	
		float dt = (clock() - timer0)/((double)CLOCKS_PER_SEC)*1000;
		printf("Took %.2f ms, obj: %f, duration: %f \n",dt,VAL(obj_sum),(tf-t0)/60./60.);
	}
}

void saveSpaceshot() {
	std::default_random_engine gen;
	gen.seed(1);
	std::normal_distribution<float> init_altitude_N(5966,10);
	std::normal_distribution<float> decent_rate_N(-22.5,5);
	std::normal_distribution<float> equil_alt_N(28700,300);
	std::normal_distribution<float> ascent_rate_N(3.7,0.05);
	std::normal_distribution<float> cutdown_time_N(2.5*60*60,60*2);
	load_data(get_data_path("/proc/gfs_pred_0deg5/20181020_18"), 1500000000,1600000000);
	int t0 = 1540060749;
	int dt_wp = 60*1;
	int tmax = 48*60*60;
	int N_wp = tmax/dt_wp;
	

	for(int j = 0; j < 10000; j++){
	float init_altitude = init_altitude_N(gen);  //m 
	float ascent_rate = ascent_rate_N(gen);		// m/s
	float decent_rate = decent_rate_N(gen); 	// m/s
	float equil_alt = equil_alt_N(gen); 		// meters
	int cutdown_time = cutdown_time_N(gen); 	// seconds
	if(equil_alt > 29000) equil_alt = 29000;
	if(decent_rate > -0.1) decent_rate = -0.1;

	int start_time = -init_altitude/ascent_rate;
	float waypoints[N_wp];
	for(int i = 0; i < N_wp; i++){
		float alt = init_altitude + ascent_rate*i*dt_wp;
		if(alt > equil_alt) alt = equil_alt;
		if((i*dt_wp - start_time) > cutdown_time){
			alt =  equil_alt + ((i*dt_wp- start_time) - cutdown_time)*decent_rate;
			if(alt<400){
				tmax = (i-1)*dt_wp;
		//		printf("flight time: %f\n",(tmax-start_time)/60./60.);
				break;
			}
		}
		//printf("%f\n",alt);
		waypoints[i] = float(alt2p(alt));
	}
	WaypointController<float> pres(t0, dt_wp, waypoints);
	Simulation<float> sim(pres, j);
	sim.tmax = tmax;
	sim.run(t0, 36.89262,-121.45095+360);
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
	stocasticGradients();
	//saveSpaceshot();
	//stocasticGradients();
}
