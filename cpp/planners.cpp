#include "planners.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "Utils.h"
StochasticMPC::StochasticMPC(const char* input_db, sim_state<float> state0)  
	//: objfn(0, 360)
{
	data.load_data(input_db, 1500000000,1600000000);
	conf.state0 = state0;
}

TemporalParameters<float>  StochasticMPC::run(){
	TemporalParameters<float> params(conf.state0.t, conf.init_cmd_dt, conf.tmax, 15000, 2000);
	for(int start = 0; start<conf.n_starts; start++){
		params.resample(conf.init_cmd_dt);
		params.rand_sets(12000,16000);
		params.randn_tols();
		params.resample(conf.cmd_dt);
		float lr_set = hparams.lr_set;
		float lr_tol = hparams.lr_tol;
		float smoothed_obj;
		float smoothed_obj_last;
		bool is_fist_run = true;
		bool converged = false;
		AdjustableLowpass obj_filt(hparams.obj_filter_corner,0.5,1);
		if (save_to_file) {
			char path[PATH_MAX];
			snprintf(path, PATH_MAX, "../ignored/sim/opt.%03d.bin", start);
			file = fopen(path, "wb");
			ensure(file != 0);
			ensure(setvbuf(file, 0, _IOFBF, 16384) == 0);
		}
		for (int it=0; it<conf.n_iters_max; it++){
			TIMEIT("Gradient Step",false,
			ctrl_cmd<double>* gradients = new ctrl_cmd<double>[conf.n_samples*params.N];
			sim_state<double>* final_states = new sim_state<double>[conf.n_samples];
			double* obj_vals = new double[conf.n_samples];
			Scheduler<float> sched(-1, conf.n_samples);
			for (int run=0; run<conf.n_samples; run++) {
				sched.add([&,run]() {
					ctrl_cmd<double>* run_gradient = gradients+run*params.N;
					adept::Stack stack;
					TemporalParameters<adouble> run_params(conf.state0.t, conf.cmd_dt, conf.tmax, 12000, 2000);
					run_params = params;
					stack.new_recording();
					StochasticControllerApprox<adouble> controller(run_params, rand());
					LinInterpWind<adouble> wind(data);
		            wind.sigma = 0;
					//FinalLongitude<adouble> objfn;
					//MinDistanceToPoint<adouble> objfn(15.638795, 16.970971+360);
					MinDistanceToPoint<adouble> objfn(57.332174, -4.663537+360);
					EulerIntBal<adouble> in;
					int fname = -1;
					if (conf.write_files && (it == 0 || it == conf.n_iters_max-1 || converged)) { /*printf("saving!\n");*/ fname = conf.fname_offset + start*conf.n_samples*2 + run + conf.n_samples*(it == 0); }
					Simulation<adouble> sim(controller, wind, objfn, in, fname);
					sim.tmax=conf.tmax;
					sim_state<adouble> sf = conf.state0.cast<adouble>();
					sim.run(sf);
					adouble obj_val = conf.opt_sign*objfn.getObjective();
					
					obj_val.set_gradient(1.0);
					obj_vals[run] = VAL(obj_val);
					final_states[run] = sf.cast<double>();
					stack.compute_adjoint();
					for(int j = 0; j<params.N; j++){
						run_gradient[j].tol = run_params.cmds[j].tol.get_gradient(); 
						run_gradient[j].h = run_params.cmds[j].h.get_gradient();
					}
					//printf("run:%d,p:%p\n",run,run_gradient);
					//printf("f;%f,%p,%d\n",(obj_vals[run]),&(obj_vals[run]),run);
					return 0.;
				});
			}
			sched.finish();
			if(converged) break; 
			ctrl_cmd<double> gradient[params.N];
			for(int j = 0; j<params.N; j++){
				gradient[j].h = 0;
				gradient[j].tol = 0;
			}
			float meanobj=0;
			float meanbal=0;
			float meantime=0;
			for (int run=0; run<conf.n_samples; run++){
				for(int j = 0; j<params.N; j++){
					gradient[j].h -= gradients[run*params.N+j].h/conf.n_samples;
					gradient[j].tol -= gradients[run*params.N+j].tol/conf.n_samples;
					//printf("d[%d],%d = (%f,%f)\n",j,run,gradients[run*params.N+j].h,gradients[run*params.N+j].tol);
				}
				meanobj += obj_vals[run]/conf.n_samples;
				meanbal += final_states[run].bal/conf.n_samples;
				meantime += (final_states[run].t-conf.state0.t)/(float)conf.n_samples;		
			}

			lr_set = hparams.alpha*lr_set;
			lr_tol = hparams.alpha*lr_tol;
			if(!is_fist_run){
				smoothed_obj = obj_filt.update(meanobj);
				if(smoothed_obj - smoothed_obj_last > 0 && it >= conf.n_iters_min) {
					converged = true;
					printf("Converged after %d iterations. obj: %f\n",it,meanobj);
				}
				smoothed_obj_last = smoothed_obj;
			} else {
				smoothed_obj = meanobj;
				obj_filt.setSS(meanobj);
				is_fist_run = false;

			}
			for(int j = 0; j<params.N; j++){
			//	printf("d[%d] = (%f,%f)\n",j,gradient[j].h,gradient[j].tol);
				double set = params.cmds[j].h + lr_set * gradient[j].h;
				set = min(hparams.set_max, max(hparams.set_min, set));
				params.cmds[j].h = set;

				double tol = params.cmds[j].tol + lr_tol * gradient[j].tol;
				tol = min(hparams.tol_max, max(hparams.tol_min, tol));
				params.cmds[j].tol = tol;
			}
			if (save_to_file) {
				fwrite(&meanobj, sizeof(float), 1, file);
				fwrite(&lr_set, sizeof(float), 1, file);
				fwrite(&smoothed_obj, sizeof(float), 1, file);
			}
			)
			//printf("obj: %f, bal: %f, days: %f\n",meanobj, meanbal, meantime/86400.);
			delete gradients;
			delete final_states;
			delete obj_vals;
		}
	}
	//printf("%d\n",params.T);

	if (save_to_file) fclose(file);
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
