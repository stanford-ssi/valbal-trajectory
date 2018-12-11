#include "planners.h"

StocasticMPC::StocasticMPC(const char* input_db, sim_state<float> state0)  
	//: objfn(0, 360)
{
	data.load_data(input_db, 1500000000,1600000000);
	conf.state0 = state0;
}

TemporalParameters<adouble>  StocasticMPC::run(){
	TemporalParameters<adouble> params(conf.state0.t, conf.cmd_dt, conf.state0.t+conf.tmax, 14000, 2000);
	for (int it=0; it<conf.n_iters; it++){
		clock_t timer0 = clock();
		adouble obj_sum = 0;
		stack.new_recording();
		adouble objectives[conf.n_samples];
		float meanbal = 0;
		float meantime = 0;
		for (int run=0; run<conf.n_samples; run++) {
			objfn.clear();
			StochasticControllerApprox<adouble> controller(params, rand());
			LinInterpWind<adouble> wind(data);
            wind.sigma = 0;
			//FinalLongitude<adouble> obj;
			//MinDistanceToPoint<adouble> obj(13.589181, -85.584796+360);
			EulerIntBal<adouble> in;
			int fname = -1;
			if (it == 0 || it == conf.n_iters-1) { printf("saving!\n"); fname = conf.fname_offset + run + conf.n_samples*(it == 0); }
			Simulation<adouble> sim(controller, wind, objfn, in, fname);
			sim.tmax=60*60*120;
			sim_state<adouble> sf = sim.run(conf.state0.t, conf.state0.lat, conf.state0.lon);
			meanbal += VAL(sf.bal);
			meantime += (sf.t - conf.state0.t);
			objectives[run] = conf.opt_sign*objfn.getObjective();
		}
		meanbal /= conf.n_samples;
		meantime /= conf.n_samples;

		for (int run=0; run<conf.n_samples; run++) obj_sum += objectives[run];
		obj_sum = obj_sum/((float)conf.n_samples);
		obj_sum.set_gradient(1.0);
		stack.compute_adjoint();
		params.apply_gradients(step);

		float dt = (clock() - timer0)/((double)CLOCKS_PER_SEC)*1000;
		printf("Took %.2f ms, obj: %f, bal: %f, days: %f\n",dt,VAL(obj_sum), meanbal, meantime/86400.);
	}
	//TemporalParameters<float> fparams(conf.state0.t, conf.cmd_dt, conf.state0.t+conf.tmax, 14000, 2000);
	//fparams = params
	return params;
}
