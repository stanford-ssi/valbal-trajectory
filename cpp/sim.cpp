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

template<class Float>
WaypointController<Float>::WaypointController(int t0_, int dt_, Float *alts_)
		: t0(t0_), dt(dt_), alts(alts_) {}

template<class Float>
Float WaypointController<Float>::get_pressure(int t, float lat, float lon) {
	unsigned int idx = (t-t0)/dt;
	float theta = (t - (t0 + dt*idx))/((float)dt);
	return alts[idx] + theta * (alts[idx+1] - alts[idx]);
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
Float PressureTable<Float>::get_pressure(int t, float lat, float lon) {
	unsigned int idx = (t-t0)/dt;
	assert(idx >= 0 && idx < n);
	return alts[idx];
}

template<class Float>
Float LasSim<Float>::get_pressure(int t, float lat, float lon){
	if(is_first_run){
		is_first_run = false;
	} else {
		const int N = int((t-t_last)*freq);
		sim.conf.lat = lat;
		sim.conf.lon = lon;
		LasagnaController::Input input;
		for(int i = 0; i < N; i++){
			input.h_abs = sim.evolve(double(las.getAction()));
			input.h_rel = input.h_abs;
			input.dldt_ext = sim.sunset_dldt*3;
			las.update(input);
		}
	}
	t_last = t;
	return alt2p(sim.h);
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
wind_vector<Float> Simulation<Float>::get_wind(int t, Float lat, Float lon, Float pres) {
	while (files[cur_file+1].time < t) {
		cur_file++;
	}
	assert(cur_file < num_files && "file requested out of range");
	debugf("%d cur file %d %d %d pres %f\n", t, cur_file, ((int)(t-files[cur_file].time)), (int)(files[cur_file].time), VAL(pres));

	/* Get pressure level. Here a simple linear search is faster, although it
	 * can probably be vectorized. */
	int i;
	for (i=0; i<NUM_LEVELS; i++) {
		if (pres <= LEVELS[i]) break;
	}
	assert(i > 0 && i < NUM_LEVELS);
	Float theta_pr = (pres-LEVELS[i-1])/(LEVELS[i] - LEVELS[i-1]);

	#define INTERP_ALT(dst, src, idx) \
		Float dst = src->data[i-1][idx] + theta_pr*(src->data[i][idx] - src->data[i-1][idx]);
	#define LAT(x) (LAT_MIN + LAT_D * (x))
	#define LON(x) (LON_MIN + LON_D * (x))

	point pt = get_base_neighbor(VAL(lat), VAL(lon));
	Float theta_lat = (lat - LAT(pt.lat))/LAT_D;
	Float theta_lon = (lon - LON(pt.lon))/LON_D;
	float theta_t = t-files[cur_file].time;
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
		debugf("v11! @ %f %d %d %f %f\n", VAL(pres), pt.lat, pt.lon, LAT(pt.lat), LON(pt.lon));
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
Simulation<Float>::Simulation(PressureSource<Float>& s, int i)
		: pressure(s), random_gen((std::random_device())()), normal(0,1) {
	save_to_file = i >= 0;
	if (save_to_file) {
		char path[PATH_MAX];
		snprintf(path, PATH_MAX, "../ignored/sim/output.%03d.bin", i);
		file = fopen(path, "wb");
		ensure(file != 0);
		ensure(setvbuf(file, 0, _IOFBF, 16384) == 0);
	}
}

template<class Float>
vec2<Float> Simulation<Float>::run(int t, Float lat, Float lon) {
	int Tmax = tmax + t;
	const float idlat = dt / (2 * M_PI * 6371008 / 360.);
	debugf("Starting from (%f, %f)\n", VAL(lat), VAL(lon));
	while (t < Tmax) {
		Float p = pressure.get_pressure(t, VAL(lat), VAL(lon));
		if (save_to_file) {
			float actual_lat = VAL(lat);
			float actual_lon = VAL(lon);
			fwrite(&actual_lat, sizeof(float), 1, file);
			fwrite(&actual_lon, sizeof(float), 1, file);

			float actual_p = VAL(p);

			/* Computing pressure -> altitude makes the simulation ~16% slower, so
			 * it's optional. */
			#ifdef STORE_ALTITUDE
				actual_p = p2alt(actual_p);
			#endif

			fwrite(&actual_p, sizeof(float), 1, file);
		}
		wind_vector<Float> w = get_wind(t, lat, lon, p);
		lat += w.v * idlat;
		lon += w.u * idlat / fastcos(lat * M_PI / 180.);
		t += dt;
	}
	debugf("Ended up at (%f, %f)\n", VAL(lat), VAL(lon));
	vec2<Float> ret = {lat, lon};
	if (save_to_file) {
		fclose(file);
	}
	return ret;
}

float p2alt(float p){
	return (1.0-(pow((p/101350.0),0.190284)))*145366.45*0.3048;
}

float alt2p(float alt){
	return pow(-((alt/145366.45/0.3048)-1.0),1.0/0.190284)*101350.0;
}

#define INIT(type) \
		template class PressureTable<type>; \
		template class WaypointController<type>; \
		template class LasSim<type>; \
		template class Simulation<type>;

#include <adept.h>
INIT(adept::adouble)
INIT(float)

