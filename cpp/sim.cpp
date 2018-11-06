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
	if(is_first_run){
		is_first_run = false;
	} else {
		const int N = int((state.t-t_last)*freq);
		sim.conf.lat = VAL(state.lat);
		sim.conf.lon = VAL(state.lon);
		LasagnaController::Input input;
		for(int i = 0; i < N; i++){
			input.h_abs = sim.evolve(double(las.getAction()));
			input.h_rel = input.h_abs;
			input.dldt_ext = sim.sunset_dldt*3;
			las.update(input);
		}
	}
	t_last = state.t;
}

template<class Float>
LasSim<Float>::LasSim(int seed, float h, float l) : las(this->freq), sim(seed) {
	sim.h = h;
	sim.l = l;
	sim.conf.freq = freq;
}

template<class Float>
LasSim<Float>::LasSim(int seed,float h) : las(this->freq), sim(seed) {
	sim.h = h;
	sim.conf.freq = freq;
}

template<class Float>
LasSim<Float>::LasSim(int seed) : las(this->freq), sim(seed) {
	sim.conf.freq = freq;
}

template<class Float>
StocasticControllerApprox<Float>::StocasticControllerApprox(int t0_, int dt_, ctrl_cmd<Float> *cmds_, int seed) : 
las_sim(seed), t0(t0_), dt(dt_), cmds(cmds_) 
{ 
	h_mid = las_sim.las.getConstants().h_cmd;
	tol0 = las_sim.las.getConstants().ss_error_thresh;
}

template<class Float>
void StocasticControllerApprox<Float>::get_pressure(sim_state<Float>& state){
	unsigned int idx = (state.t-t0)/dt;
	float theta = (state.t - (t0 + dt*idx))/((float)dt);
	Float cmd_h = cmds[idx].h + theta * (cmds[idx+1].h - cmds[idx].h);
	Float cmd_tol = cmds[idx].tol + theta * (cmds[idx+1].tol - cmds[idx].tol);
	state.bal_rate = 0.03/60./60.*750/cmd_tol + 0.03/60./60.;  
	sim_state<float> state_val;
	state_val.lat = VAL(state.lat);
	state_val.lon = VAL(state.lon);
	state_val.t = state.t;
	las_sim.get_altitude(state_val);
	float las_h = state_val.p;
	Float alt = cmd_h + (las_h - h_mid)*cmd_tol/tol0;
	state.p = alt2p(alt);
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
	while (files[cur_file+1].time < s.t) {
		cur_file++;
	}
	assert(cur_file < num_files && "file requested out of range");
	debugf("%d cur file %d %d %d pres %f\n", s.t, cur_file, ((int)(s.t-files[cur_file].time)), (int)(files[cur_file].time), VAL(s.p));

	/* Get pressure level. Here a simple linear search is faster, although it
	 * can probably be vectorized. */
	int i;
	for (i=0; i<NUM_LEVELS; i++) {
		if (s.p <= LEVELS[i]) break;
	}
	assert(i > 0 && i < NUM_LEVELS);
	Float theta_pr = (s.p-LEVELS[i-1])/(LEVELS[i] - LEVELS[i-1]);

	#define INTERP_ALT(dst, src, idx) \
		Float dst = src->data[i-1][idx] + theta_pr*(src->data[i][idx] - src->data[i-1][idx]);
	#define LAT(x) (LAT_MIN + LAT_D * (x))
	#define LON(x) (LON_MIN + LON_D * (x))

	point pt = get_base_neighbor(VAL(s.lat), VAL(s.lon));
	Float theta_lat = (s.lat - LAT(pt.lat))/LAT_D;
	Float theta_lon = (s.lon - LON(pt.lon))/LON_D;
	float theta_t = s.t-files[cur_file].time;
	theta_t /= (files[cur_file+1].time-files[cur_file].time);
	debugf("theta t %f\n", theta_t);
	assert(theta_t >= 0 && theta_t <= 1);

	data_file *f = files + cur_file;
	Float us[2];
	Float vs[2];
	#ifdef INTERPOLATE_IN_TIME
	for (int j=0; j<2; j++) {
	#else
	{ int j = 0;
	#endif
		wind_data *p11 = get_data_at_point(f+j, {pt.lat, pt.lon});
		wind_data *p12 = get_data_at_point(f+j, {pt.lat, pt.lon + 1});
		wind_data *p21 = get_data_at_point(f+j, {pt.lat + 1, pt.lon});
		wind_data *p22 = get_data_at_point(f+j, {pt.lat + 1, pt.lon + 1});

		INTERP_ALT(u11, p11, 0);
		INTERP_ALT(u21, p21, 0);
		Float ulat1 = u11 + theta_lat * (u21 - u11);
		debugf("v11! @ %f %d %d %f %f\n", VAL(s.p), pt.lat, pt.lon, LAT(pt.lat), LON(pt.lon));
		for (int k=0; k<NUM_LEVELS; k++) {
			debugf("%f: %d\n", LEVELS[k], p11->data[k][1]);
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
		: pressure(s), wind(w), intg(in), objfn(o) {
	init(i);
}

template<class Float>
Simulation<Float>::Simulation(PressureSource<Float>& s, WindSource<Float>& w, ObjectiveFn<Float>& o, int i)
		: intg_default(), pressure(s), wind(w), intg(intg_default), objfn(o) {
	init(i);
}

template<class Float>
Simulation<Float>::Simulation(PressureSource<Float>& s, int i)
		: wind_default(), obj_default(), intg_default(), pressure(s), wind(wind_default), intg(intg_default), objfn(obj_default) {
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
sim_state<Float> Simulation<Float>::run(int t, Float lat, Float lon) {
	int Tmax = tmax + t;
	sim_state<Float> state;
	state.lat = lat;
	state.lon = lon;
	state.t = t;
	state.bal = 10;
	debugf("Starting from (%f, %f)\n", VAL(lat), VAL(lon));
	while (state.t < Tmax) {
		pressure.get_pressure(state);
		if (save_to_file) {
			float actual_lat = VAL(state.lat);
			float actual_lon = VAL(state.lon);
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
	return state;
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
		template class StocasticControllerApprox<type>;

#include <adept.h>
INIT_SIM(adept::adouble)
INIT_SIM(float)

