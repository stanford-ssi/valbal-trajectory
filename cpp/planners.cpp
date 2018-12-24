#include "planners.h"

StochasticMPC::StochasticMPC(const char* input_db, sim_state<float> state0)  
	//: objfn(0, 360)
{
	data.load_data(input_db, 1500000000,1600000000);
	conf.state0 = state0;
}

TemporalParameters<float>  StochasticMPC::run(){
	TemporalParameters<float> params(conf.state0.t, conf.cmd_dt, conf.tmax, 12000, 2000);
	for (int it=0; it<conf.n_iters; it++){
		TIMEIT("Gradient Step",
		ctrl_cmd<double> gradients[conf.n_samples][params.N];
		sim_state<double> final_states[conf.n_samples];
		double obj_vals[conf.n_samples];
		Scheduler<float> sched(-1, conf.n_samples);
		for (int run=0; run<conf.n_samples; run++) {
			//ctrl_cmd<double>* run_gradient = &gradients[run][0];
			sched.add([&,run]() {
				/*
				adept::Stack stack;
				TemporalParameters<adouble> run_params(conf.state0.t, conf.cmd_dt, conf.tmax, 12000, 2000);
				run_params = params;
				stack.new_recording();
				StochasticControllerApprox<adouble> controller(run_params, rand());
				LinInterpWind<adouble> wind(data);
	            wind.sigma = 0;
				//FinalLongitude<adouble> objfn;
				MinDistanceToPoint<adouble> objfn(15.638795, 16.970971+360);
				EulerIntBal<adouble> in;
				int fname = -1;
				if (conf.write_files && (it == 0 || it == conf.n_iters-1)) { printf("saving!\n"); fname = conf.fname_offset + run + conf.n_samples*(it == 0); }
				Simulation<adouble> sim(controller, wind, objfn, in, fname);
				sim.tmax=conf.tmax;
				sim_state<adouble> sf = conf.state0.cast<adouble>();
				sim.run(sf);
				adouble obj_val = conf.opt_sign*objfn.getObjective();
				*/
				//obj_val.set_gradient(1.0);
				obj_vals[run] = run*0.34;//VAL(obj_val);
				//final_states[run] = sf.cast<double>();
				//stack.compute_adjoint();
				//for(int j = 0; j<params.N; j++){
				//	run_gradient[j].tol = run_params.cmds[j].tol.get_gradient(); 
				//	run_gradient[j].h = run_params.cmds[j].h.get_gradient();
				//}
				printf("f;%f,%p,%d\n",(obj_vals[run]),&(obj_vals[run]),run);
				return 0.;
			});
		}
		sched.finish();
		ctrl_cmd<double> gradient[params.N];
		double meanobj=0;
		double meanbal=0;
		double meantime=0;
		for (int run=0; run<conf.n_samples; run++){
			for(int j = 0; j<params.N; j++){
				gradient[run].h += gradients[run][j].h/conf.n_samples;
				gradient[run].tol += gradients[run][j].tol/conf.n_samples;
				//printf("d(%f,%f)\n",gradient[run].h,gradient[run].tol);
			}
			meanobj += obj_vals[run]/conf.n_samples;
			printf("g;%f,%p,%d\n",(obj_vals[run]),&(obj_vals[run]),run);
			meanbal += final_states[run].bal/conf.n_samples;
			meantime += (final_states[run].t-conf.state0.t)/(float)conf.n_samples;
		}
		)
		printf("obj: %f, bal: %f, days: %f\n",meanobj, meanbal, meantime/86400.);
	}
	//printf("%d\n",params.T);
	return params;
}


SpatialPlanner::SpatialPlanner(const char* input_db, sim_state<float> state0)  
	//: objfn(0, 360)
{
	data.load_data(input_db, 1500000000,1600000000);
	conf.state0 = state0;
}

SpatiotemporalParameters<float>  SpatialPlanner::run(){
	SpatiotemporalParameters<adouble> params(conf.state0.t, conf.cmd_dt, conf.tmax + conf.cmd_dt, 12000, 2000);
	for (int it=0; it<conf.n_iters; it++){
		clock_t timer0 = clock();
		adouble obj_sum = 0;
		stack.new_recording();
		adouble objectives[conf.n_samples];
		float meanbal = 0;
		float meantime = 0;
		for (int run=0; run<conf.n_samples; run++) {
			//MinDistanceToPoint<adouble> objfn(46.225336, -74.891043+360);
			StochasticControllerApprox<adouble> controller(params, rand());
			LinInterpWind<adouble> wind(data);
            wind.sigma = 0;
			FinalLongitude<adouble> objfn;
			//MinDistanceToPoint<adouble> obj(13.589181, -85.584796+360);
			EulerIntBal<adouble> in;
			int fname = -1;
			if (conf.write_files && (it == 0 || it == conf.n_iters-1)) { printf("saving!\n"); fname = conf.fname_offset + run + conf.n_samples*(it == 0); }
			Simulation<adouble> sim(controller, wind, objfn, in, fname);
			sim.tmax=conf.tmax;
			sim_state<adouble> sf = conf.state0.cast<adouble>();
			sim.run(sf);
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
	SpatiotemporalParameters<float> fparams(conf.state0.t, conf.cmd_dt, conf.tmax, 14000, 2000);
	fparams = params;
	//printf("%d\n",params.T);
	return fparams;
}
