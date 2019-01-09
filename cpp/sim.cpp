#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <random>

#include "sim.h"
#include "utils.h"

#define STORE_ALTITUDE

template <class Float>
void EulerInt<Float>::integrate(sim_state<Float>& loc, wind_vector<Float>& w){
	loc.lat += w.v * idlat;
	loc.lon += w.u * idlat / fastcos(loc.lat * M_PI / 180.);
	loc.t += dt;
}

template <class Float>
void EulerIntBal<Float>::integrate(sim_state<Float>& s, wind_vector<Float>& w){
	if(s.bal <= 0){
		s.t += dt;
		return;
	}
	Float bal = s.bal - s.bal_rate*(float)dt ;
	if(bal <= 0){
		Float dts = s.bal/s.bal_rate;
		s.lat += w.v * dts/(float)dt * idlat;
		s.lon += w.u * dts/(float)dt * idlat / fastcos(s.lat * M_PI / 180.);
	} else {
		s.lat += w.v * idlat;
		s.lon += w.u * idlat / fastcos(s.lat * M_PI / 180.);
	}
	s.bal = bal;
	s.t += dt;

}

template<class Float>
WaypointController<Float>::WaypointController(int t0_, int dt_, Float *alts_)
		: t0(t0_), dt(dt_), alts(alts_) {}

template<class Float>
void WaypointController<Float>::get_pressure(sim_state<Float>& state) {
	unsigned int idx = (state.t-t0)/dt;
	float theta = (state.t - (t0 + dt*idx))/((float)dt);
	state.p = alts[idx] + theta * (alts[idx+1] - alts[idx]);
}

template<class Float>
PressureTable<Float>::PressureTable(const char *fname) {
	FILE *f = fopen(fname, "rb");
	assert(f != 0);

	ensure(fread(&dt, 4, 1, f) == 1);
	ensure(fread(&t0, 4, 1, f) == 1);
	ensure(fread(&n, 4, 1, f) == 1);

	alts = (float*)malloc(sizeof(float)*n);
	ensure(fread(alts, sizeof(float), n, f) == n);
	ensure(fclose(f) == 0);

	debugf("Loaded %d pressures from %s, dt: %d s, t0: %d.\n", n, fname, dt, t0);
}

template<class Float>
void PressureTable<Float>::get_pressure(sim_state<Float>& state) {
	unsigned int idx = (state.t-t0)/dt;
	assert(idx >= 0 && idx < n);
	state.p = alts[idx];
}

template<class Float>
void LasSim<Float>::get_pressure(sim_state<Float>& state){
	evolve(state);
	state.p = alt2p(sim.h);
}

template<class Float>
void LasSim<Float>::get_altitude(sim_state<Float>& state){
	evolve(state);
	state.p = sim.h;
}

template<class Float>
void LasSim<Float>::evolve(sim_state<Float>& state){
	if(!(t_last || dt)){
	} else {
		const int N = dt ? int(dt*freq) : int((state.t-t_last)*freq);
		sim.conf.lat = VAL(state.lat);
		sim.conf.lon = VAL(state.lon);
		LasagnaController::Input input;
		for(int i = 0; i < N; i++){
			if(changing_cmds){
				sim_state<float> statef = state.template cast<float>();
				lasconst.setpoint = cmds.get_param(statef).h;
				lasconst.tolerance = cmds.get_param(statef).tol;
				las.updateConstants(lasconst);
				state.bal_rate = 0.03/60./60.*750/lasconst.tolerance + 0.04/60./60.;
				//printf("[cmds] t:%f tol:%f set:%f\n",(state.t-1543492801)/60./60.,lasconst.tolerance,lasconst.setpoint);
			}
			input.h_abs = sim.evolve(double(las.getAction()));
			input.h_rel = input.h_abs;
			input.dldt_ext = sim.sunset_dldt*3;
			las.update(input);
		}
		state.cmd.h = lasconst.setpoint;
		state.cmd.tol = lasconst.tolerance;
	}
	t_last = state.t;
}

template<class Float>
LasSim<Float>::LasSim(int seed, float h, float l,TemporalParameters<float>& cmds) : las(this->freq), sim(seed), cmds_defualt(0,1,0,1.,1.), cmds(cmds){
	sim.h = h;
	sim.l = l;
	sim.conf.freq = freq;
	changing_cmds = true;
}

template<class Float>
LasSim<Float>::LasSim(int seed, float h, float l) : las(this->freq), sim(seed), cmds_defualt(0,1,0,1.,1.), cmds(cmds_defualt) {
	sim.h = h;
	sim.l = l;
	sim.conf.freq = freq;
}


template<class Float>
LasSim<Float>::LasSim(int seed,float h) : las(this->freq), sim(seed), cmds_defualt(0,1,0,1.,1.), cmds(cmds_defualt) {
	sim.h = h;
	sim.conf.freq = freq;
}

template<class Float>
LasSim<Float>::LasSim(int seed) : las(this->freq), sim(seed),cmds_defualt(0,1,0,1.,1.), cmds(cmds_defualt)  {
	sim.conf.freq = freq;
}

template<class Float>
StochasticControllerApprox<Float>::StochasticControllerApprox(ParameterServer<Float>& ps_, int seed) : 
las_sim(seed), params(ps_)
{ 
	h_mid = las_sim.las.getConstants().setpoint;
	tol0 = las_sim.las.getConstants().tolerance;
	//TODO: Load ballast coeffs from file
}

template<class Float>
Float StochasticControllerApprox<Float>::get_ballast_rate(Float tol){
	tol = tol/1000; // convert to km (should eventually have everything be km)
	Float bal = 0;
	for(int i = 0;i<8;i++){
		bal += bal_coeffs[i]*pow(tol,7-i);
	}
	return bal;
}

template<class Float>
void StochasticControllerApprox<Float>::get_pressure(sim_state<Float>& state){
	ctrl_cmd<Float> cmd = params.get_param(state);
	Float cmd_h = cmd.h;
	Float cmd_tol = cmd.tol;
	state.bal_rate = get_ballast_rate(cmd_tol); 
	
	//TODO: Add nighfall ballast use

	sim_state<float> state_val;
	state_val.lat = VAL(state.lat);
	state_val.lon = VAL(state.lon);
	state_val.t = state.t;
	las_sim.get_altitude(state_val);
	float las_h = state_val.p;
	Float alt;
	if(fabs(las_h - h_mid) > tol0){
		alt = cmd_h + sgn(las_h - h_mid)*(cmd_tol + (fabs(las_h - h_mid) - tol0));
	} else {
		alt = cmd_h + (las_h - h_mid)*cmd_tol/tol0;
	}
	state.cmd = cmd; 
	state.p = alt2p(alt);
}

template<class Float>
TemporalParameters<Float>::TemporalParameters(int t0_, int dt_, int T_, double default_h_, double default_tol_) : 
t0(t0_), dt(dt_), T(T_), default_h(default_h_), default_tol(default_tol_)
{ 
	N = ceil(T_/dt_) + 1; 
	cmds = new ctrl_cmd<Float>[N];
	for (int i=0; i<N; i++) {
		cmds[i].h = default_h_;
		cmds[i].tol = default_tol_;
	}
}

template<class Float>
TemporalParameters<Float>::~TemporalParameters() {
	delete[] cmds;
}

template<class Float>
void TemporalParameters<Float>::rand_sets(double min_, double max_){
    std::mt19937 random_gen(time(0));
    std::uniform_real_distribution<double> dist(min_,max_);
    for (int i=0; i<N; i++) {
    	cmds[i].h = dist(random_gen);
    };
}

template<class Float>
void TemporalParameters<Float>::rand_tols(double min_, double max_){
    std::mt19937 random_gen(time(0));
    std::uniform_real_distribution<double> dist(min_,max_);
    for (int i=0; i<N; i++) {
    	cmds[i].tol = dist(random_gen);
    };
}

template<class Float>
void TemporalParameters<Float>::randn_tols(double mean_, double std_, double min_, double max_){
    std::mt19937 random_gen(time(0));
    std::normal_distribution<double> dist(mean_,std_);
    for (int i=0; i<N; i++) {
    	double val = dist(random_gen);
    	cmds[i].tol = min(max_, max(min_, val));
    };
}

template<class Float>
void TemporalParameters<Float>::resample(int new_dt){
	int new_N = ceil(T/new_dt) + 1;
	ctrl_cmd<Float>* new_cmds = new ctrl_cmd<Float>[new_N];
	sim_state<Float> state;
	state.t = t0;
	for (int i=0; i<new_N; i++) {
		new_cmds[i] = this->get_param(state);
		state.t += new_dt;
		//printf("set:%f\n",VAL(new_cmds[i].h));
	}
	delete[] cmds;
	cmds = new_cmds;
	dt = new_dt;
	N = new_N;
}


template<class Float>
ctrl_cmd<Float> TemporalParameters<Float>::get_param(sim_state<Float>& state){
	unsigned int idx = (state.t-t0)/dt;
	float theta = (state.t - (t0 + dt*idx))/((float)dt);
	//printf("idx:%d,theta:%f,N:%d\n",idx,theta,N);
	if((int)idx+1 > N-1){theta = 1; idx = N-2;} 
	ctrl_cmd<Float> cmd;
	cmd.h = cmds[idx].h + theta * (cmds[idx+1].h - cmds[idx].h);
	cmd.tol = cmds[idx].tol + theta * (cmds[idx+1].tol - cmds[idx].tol);
	return cmd;
}

template <class Float>
double TemporalParameters<Float>::apply_gradients(StepRule& opt) {
	return apply_gradients(opt, tag<TemporalParameters>());
}

template <class Float>
double TemporalParameters<Float>::apply_gradients(StepRule& opt, tag<TemporalParameters<float>>) {
	printf("You what mate, what are you trying to take the gradient of\n");
	exit(1);
	return M_PI;
}

template <class Float>
double TemporalParameters<Float>::apply_gradients(StepRule& step_rule, tag<TemporalParameters<adouble>>) {
	ctrl_cmd<adouble> *cmds_ = (ctrl_cmd<adouble>*)(&cmds[0]);
	step_rule.new_step();
	step_rule.optimize(cmds_,(int)(T/dt));
	return 0.0;
}

template<class Float>
void GreedySearch<Float>::get_pressure(sim_state<Float>& state){
	float min=1000000;
	Float argmin = (range.b + range.a)/2;
	for(int i=0; i < N_levels; i++){
		state.p = range.a + (range.b - range.a)*(i+0.5)/N_levels;
		wind_vector<Float> w = wind.get_wind(state);
		intg.integrate(state,w);
		Float cost = objfn.update(state,false);
		if(VAL(cost) < min){
			argmin = state.p;
			min = VAL(cost);
		}
	}
	state.p = argmin;
}

template<class Float>
wind_vector<Float> LinInterpWind<Float>::get_wind(sim_state<Float>& s) {
	while (dat.files[cur_file+1].time < s.t) {
		cur_file++;
	}
	assert(cur_file < dat.num_files && "file requested out of range");
	debugf("%d cur file %d %d %d pres %f\n", s.t, cur_file, ((int)(s.t-dat.files[cur_file].time)), (int)(dat.files[cur_file].time), VAL(s.p));

	/* Get pressure level. Here a simple linear search is faster, although it
	 * can probably be vectorized. */
	int i;
	for (i=0; i<dat.NUM_LEVELS; i++) {
		if (s.p <= dat.LEVELS[i]) break;
	}
	assert(i > 0 && i < dat.NUM_LEVELS);
	Float theta_pr = (s.p-dat.LEVELS[i-1])/(dat.LEVELS[i] - dat.LEVELS[i-1]);

	#define INTERP_ALT(dst, src, idx) \
		Float dst = src[dat.NUM_VARIABLES*(i-1) + idx] + theta_pr*(src[dat.NUM_VARIABLES*i + idx] - src[dat.NUM_VARIABLES*(i-1) + idx]);
	#define LAT(x) (dat.LAT_MIN + dat.LAT_D * (x))
	#define LON(x) (dat.LON_MIN + dat.LON_D * (x))

	point pt = dat.get_base_neighbor(VAL(s.lat), VAL(s.lon));
	Float theta_lat = (s.lat - LAT(pt.lat))/dat.LAT_D;
	Float theta_lon = (s.lon - LON(pt.lon))/dat.LON_D;
	float theta_t = s.t-dat.files[cur_file].time;
	theta_t /= (dat.files[cur_file+1].time-dat.files[cur_file].time);
	debugf("theta t %f\n", theta_t);
	assert(theta_t >= 0 && theta_t <= 1);

	data_file *f = dat.files + cur_file;
	Float us[2];
	Float vs[2];
	#ifdef INTERPOLATE_IN_TIME
	for (int j=0; j<2; j++) {
	#else
	{ int j = 0;
	#endif

		wind_t *p11 = dat.get_data_at_point(f+j, {pt.lat, pt.lon});
		wind_t *p12 = dat.get_data_at_point(f+j, {pt.lat, pt.lon + 1});
		wind_t *p21 = dat.get_data_at_point(f+j, {pt.lat + 1, pt.lon});
		wind_t *p22 = dat.get_data_at_point(f+j, {pt.lat + 1, pt.lon + 1});

		INTERP_ALT(u11, p11, 0);
		INTERP_ALT(u21, p21, 0);
		Float ulat1 = u11 + theta_lat * (u21 - u11);
		debugf("v11! @ %f %d %d %f %f\n", VAL(s.p), pt.lat, pt.lon, LAT(pt.lat), LON(pt.lon));
		for (int k=0; k<dat.NUM_LEVELS; k++) {
			debugf("%f: %d\n", dat.LEVELS[k], p11[dat.NUM_LEVELS*k+1]);
		}

		INTERP_ALT(u12, p12, 0);
		INTERP_ALT(u22, p22, 0);
		Float ulat2 = u12 + theta_lat * (u22 - u12);

		us[j] = ulat1 + theta_lon * (ulat2 - ulat1);
		debugf("U: %f (%f,%f,%f,%f)\n", VAL(us[j]),
				VAL(u11), VAL(u12), VAL(u21), VAL(u22));

		INTERP_ALT(v11, p11, 1);
		INTERP_ALT(v21, p21, 1);
		Float vlat1 = v11 + theta_lat * (v21 - v11);

		INTERP_ALT(v12, p12, 1);
		INTERP_ALT(v22, p22, 1);
		Float vlat2 = v12 + theta_lat * (v22 - v12);

		vs[j] = vlat1 + theta_lon * (vlat2 - vlat1);
		debugf("V: %f (%f,%f,%f,%f)\n", VAL(vs[j]),
				VAL(v11), VAL(v12), VAL(v21), VAL(v22));
	}
	debugf("---\n");

	/* When I find myself in need of variance,
	 * Father Wareham comes to me,
	 * Speaking words of cunning,
	 * sqrt(E[X^2] - E[X]^2) = es. tee. dee */

	wind_vector<Float> w;
	#ifdef INTERPOLATE_IN_TIME
		w.u = (us[0] + theta_t * (us[1] - us[0])) / 100.;
		w.v = (vs[0] + theta_t * (vs[1] - vs[0])) / 100.;
	#else
		w.u = us[0];
		w.v = vs[0];
	#endif
	if (sigma != 0) {
		w.u += sigma * normal(random_gen);
		w.v += sigma * normal(random_gen);
	}

	return w;
}

template<class Float>
Simulation<Float>::Simulation(PressureSource<Float>& s, WindSource<Float>& w, ObjectiveFn<Float>& o, Integrator<Float>& in, int i)
		: wind_default(dat_default), pressure(s), wind(w), intg(in), objfn(o) {
	init(i);
}

template<class Float>
Simulation<Float>::Simulation(PressureSource<Float>& s, WindSource<Float>& w, ObjectiveFn<Float>& o, int i)
		: wind_default(dat_default), intg_default(), pressure(s), wind(w), intg(intg_default), objfn(o) {
	init(i);
}

template<class Float>
Simulation<Float>::Simulation(PressureSource<Float>& s, DataHandler& d, int i)
		: wind_default(d), obj_default(), intg_default(), pressure(s), wind(wind_default), intg(intg_default), objfn(obj_default) {
	init(i);
}

template<class Float>
void Simulation<Float>::init(int i){
	save_to_file = i >= 0;
	calc_obj = (dynamic_cast<NoOp<Float>*> (&objfn) == NULL);
	if (save_to_file) {
		char path[PATH_MAX];
		snprintf(path, PATH_MAX, "../ignored/sim/output.%03d.bin", i);
		file = fopen(path, "wb");
		ensure(file != 0);
		ensure(setvbuf(file, 0, _IOFBF, 16384) == 0);
	}
}

template<class Float>
void Simulation<Float>::run(sim_state<Float>& state) {
	debugf("Starting from (%f, %f)\n", VAL(state.lat), VAL(state.lon));
	int Tmax = tmax + state.t;
	while (state.t < Tmax) {
		pressure.get_pressure(state);
		if (save_to_file) {
			float actual_lat = VAL(state.lat);
			float actual_lon = VAL(state.lon);
			float actual_set = VAL(state.cmd.h);
			float actual_tol = VAL(state.cmd.tol);
			fwrite(&actual_lat, sizeof(float), 1, file);
			fwrite(&actual_lon, sizeof(float), 1, file);

			float actual_p = VAL(state.p);

			/* Computing pressure -> altitude makes the simulation ~16% slower, so
			 * it's optional. */
			#ifdef STORE_ALTITUDE
				actual_p = p2alt(actual_p);
			#endif
			debugf("[sim state] time:%d lat:%.4f lon:%.4f alt:%.1f bal:%.2f\n",state.t, actual_lat, actual_lon, actual_p, VAL(state.bal));

			fwrite(&actual_p, sizeof(float), 1, file);
			fwrite(&actual_set, sizeof(float), 1, file);
			fwrite(&actual_tol, sizeof(float), 1, file);
		}
		wind_vector<Float> w = wind.get_wind(state);
		intg.integrate(state,w);
		if(calc_obj) objfn.update(state);
		if(state.bal < 0) break;
	}
	debugf("Ended up at (%f, %f)\n", VAL(state.lat), VAL(state.lon));
	if (save_to_file) {
		fclose(file);
	}
}

template<class Float>
sim_state<Float> Simulation<Float>::run(int t, Float lat, Float lon) {
	sim_state<Float> state;
	state.lat = lat;
	state.lon = lon;
	state.t = t;
	state.bal = 4.5;
	run(state);
	return state;
}

template<class Float>
SpatiotemporalParameters<Float>::SpatiotemporalParameters(int t0_, int dt_, int T_, double default_h_, double default_tol_) : 
t0(t0_), dt(dt_), T(T_), default_h(default_h_), default_tol(default_tol_)
{ 
	int N = T_/dt_;
	cmds = new cmd_tree<Float>[N];
	for (int i=0; i<N; i++) {
		cmds[i].cmd.h = default_h_;
		cmds[i].cmd.tol = default_tol_;
	}
}

template<class Float>
SpatiotemporalParameters<Float>::~SpatiotemporalParameters() {
	delete[] cmds;
}

template<class Float>
ctrl_cmd<Float> SpatiotemporalParameters<Float>::get_param(sim_state<Float>& state){
	unsigned int idx = (state.t-t0)/dt;
	//printf("heyyy %d %f %f\n", state.t, VAL(state.lat), VAL(state.lon));
	float theta = (state.t - (t0 + dt*idx))/((float)dt);
	ctrl_cmd<Float> cmd;
	ctrl_cmd<Float>& cmd_0 = cmds[idx].get_cmd(state);
	ctrl_cmd<Float>& cmd_1 = cmds[idx+1].get_cmd(state);
	cmd.h = cmd_0.h + theta * (cmd_1.h - cmd_0.h);
	cmd.tol = cmd_0.tol + theta * (cmd_1.tol - cmd_0.tol);
	return cmd;
}

template<class Float>
ctrl_cmd<Float>& cmd_tree<Float>::get_cmd(sim_state<Float>& state) {
	if (upper == 0 || lower == 0) {
		//printf("Made it to bottom of tree! %f %f\n", VAL(state.lat), VAL(state.lon));
		array<float, D> arr;
		arr[0] = VAL(state.lat);
		arr[1] = VAL(state.lon);
		requests.push_back(arr);
		for (int i=0; i<D; i++) {
			mins[i] = min(mins[i], arr[i]);
			maxs[i] = max(maxs[i], arr[i]);
		}
		return cmd;
	}
	double v = a * VAL(state.lat) + b * VAL(state.lon);
	if (v >= c) {
		return upper->get_cmd(state);	
	} else {
		return lower->get_cmd(state);
	}
}

template<int D>
float L2dist(array<float, D> a, array<float, D> b) {
	float o = 0;
	for (int i=0; i<D; i++) {
		o += (a[i]-b[i])*(a[i]-b[i]);
	}
	return sqrt(o);
}

template<int D>
void kmeans(int k, vector<array<float, D>>& requests, array<float, D+1>& hyperplane) {
	assert(requests.size() >= (unsigned int)k);
	vector<int> cluster(requests.size());
	vector<array<float, D>> centers(k);
	vector<array<float, D>> newcenters(k);
	vector<int> assignments(requests.size());
	for (size_t i=0; i<requests.size(); i++) {
		assignments[i] = -1;
	}
	vector<int> ncenters(k);
	for (int i = 0; i < k; i++) {
		centers[i] = requests[i];
	}
	for (int it = 0; it<10; it++) {
		for (int i=0; i<k; i++) {
			for (int j=0; j<D; j++) {
				newcenters[i][j] = 0;
			}
			ncenters[i] = 0;
		}
		bool changed = false;
		for (size_t i=0; i<requests.size(); i++) {
			int bestcluster = -1;
			float bestL2 = 1e20;
			for (int j=0; j<k; j++) {
				float dst = L2dist<D>(requests[i], centers[j]);
				if (dst < bestL2) {
					bestcluster = j;
					bestL2 = dst;
				}
			}
			debugf("request %lu assigned to %d\n", i, bestcluster);
			for (int j=0; j<D; j++) {
				newcenters[bestcluster][j] += requests[i][j];
			}
			ncenters[bestcluster]++;
			if (bestcluster != assignments[i]) changed = true;
			assignments[i] = bestcluster;
		}
		debugf("Centers:\n");
		for (int i=0; i<k; i++) {
			debugf(" (");
			for (int j=0; j<D; j++) {
				centers[i][j] = newcenters[i][j]/ncenters[i];
				debugf(" %f", centers[i][j]);
			}
			debugf("\n");
		}
		if (!changed) {
			debugf("converged!\n");
			break;
		}
	}
	/* Voronoi time
		\sum (x_i - a_i)^2 > \sum (x_i - b_i)^2
		\sum (x_i^2 - 2*x_i*a_i + a_i^2) > \sum (x_i^2 - 2*x_i*b_i + b_i^2)
		lat * 2*(-a_1 + b_1) + lon * 2*(-a_2 + b_2) > b_1^2 + b_2^2 - a_1^2 - a_2^2
	*/
	if (D == 2) {
		hyperplane[0] = 2*(-centers[0][0] + centers[1][0]);
		hyperplane[1] = 2*(-centers[0][1] + centers[1][1]);
		hyperplane[2] = centers[1][0]*centers[1][0] + centers[1][1]*centers[1][1] \
						 - centers[0][0]*centers[0][0] - centers[0][1]*centers[0][1];
	} else {
		printf("not computing hyperplane!\n");
	}
	
}

template<class Float>
void cmd_tree<Float>::save_to_file(const char *fname) {
	FILE *f = fopen(fname, "wb");
	_save_to_file(f);
	fclose(f);
}

template<class Float>
void cmd_tree<Float>::_save_to_file(FILE *f) {
	if (upper == 0 || lower == 0) {
	fprintf(f, ":%lu", requests.size());
	for (size_t i=0; i<requests.size(); i++) {
		fprintf(f, " %f,%f", requests[i][0], requests[i][1]);
	}
	fprintf(f, "\n");
	}
	if (upper != 0 && lower != 0) {
		fprintf(f, ";%f %f %f\n", a, b, c);
		upper->_save_to_file(f);
		lower->_save_to_file(f);
	}
}

template<class Float>
void cmd_tree<Float>::gradients_and_split(StepRule& opt) {
	if (upper == 0 || lower == 0) {
		assert(sizeof(*this) == sizeof(cmd_tree<adept::adouble>));
		//debugf("Made it to bottom of tree while looking for gradidents!\n");
		opt.optimize((ctrl_cmd<adept::adouble>*)&cmd, 1);
		//printf("got %lu requests\n", requests.size());
		if (requests.size() == 0) return;
		// L1 distance
		if (((maxs[0]-mins[0]) + (maxs[1]-mins[1])) > 15) {
			debugf("time to split things up! %f %f\n", maxs[0]-mins[0], maxs[1]-mins[1]);
			array<float, 3> hyperplane;	
			for (size_t i=0; i<requests.size(); i++) {
				debugf("%f,%f ", requests[i][0], requests[i][1]);
			}
			kmeans<D>(2, requests, hyperplane);
			debugf("\n%f * lat + %f * lon >= %f\n", hyperplane[0], hyperplane[1], hyperplane[2]);
			a = hyperplane[0];
			b = hyperplane[1];
			c = hyperplane[2];
			upper = new cmd_tree;
			upper->cmd.h = cmd.h;
			upper->cmd.tol = cmd.tol;
			lower = new cmd_tree;
			lower->cmd.h = cmd.h;
			lower->cmd.tol = cmd.tol;
			debugf("just created %p and %p\n", upper, lower);
		} else {
			debugf("no need to split things up! %f %f\n", maxs[0]-mins[0], maxs[1]-mins[1]);
		}
		for (int i=0; i<D; i++) {
			mins[i] = 1e20;
			maxs[i] = -1e20;
		}
		requests.clear();
		return;
	}
	assert(requests.size() == 0); // We should only touch leaf nodes!
	upper->gradients_and_split(opt);	
	lower->gradients_and_split(opt);
}


template <class Float>
double SpatiotemporalParameters<Float>::apply_gradients(StepRule& opt) {
	return apply_gradients(opt, tag<SpatiotemporalParameters>());
}

template <class Float>
double SpatiotemporalParameters<Float>::apply_gradients(StepRule &opt, tag<SpatiotemporalParameters<float>>) {
	printf("You what mate, what are you trying to take the gradient of\n");
	exit(1);
	return M_PI;
}

template <class Float>
double SpatiotemporalParameters<Float>::apply_gradients(StepRule &opt, tag<SpatiotemporalParameters<adouble>>) {
	opt.new_step();
	for (int i=0; i<(T/dt); i++) {
		cmds[i].gradients_and_split(opt);
	}
    /********* NOT SURE WHAT THIS SHOULD DO WITH NEW SYNTAX TO I JUST DISABLED IT @JCREUS*******/
	//cmd_tree<adouble> *cmds_ = (cmd_tree<adouble>*)(&cmds[0]);
	//opt.optimize(cmds_,(int)(T/dt));
	return 0;
}


#define INIT_SIM(type) \
		template class PressureTable<type>; \
		template class WaypointController<type>; \
		template class LasSim<type>; \
		template class Simulation<type>; \
		template class LinInterpWind<type>; \
		template class NoOp<type>;\
		template class EulerInt<type>;\
		template class EulerIntBal<type>;\
		template class GreedySearch<type>;\
		template class FinalLongitude<type>; \
		template class StochasticControllerApprox<type>; \
		template class TemporalParameters<type>; \
		template class SpatiotemporalParameters<type>; \
		template class cmd_tree<type>;

#include <adept.h>
INIT_SIM(adept::adouble)
INIT_SIM(float)

