#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include <adept.h>
using adept::adouble;

#include "sim.h"

#define Assert(x) { if (!(x)) { printf("fatal error\n"); exit(1); } }

template<class Float>
WaypointController<Float>::WaypointController(int t0_, int dt_, Float *alts_)
		: t0(t0_), dt(dt_), alts(alts_) {}

template<class Float>
Float WaypointController<Float>::get_pressure(int t) {
	unsigned int idx = (t-t0)/dt;
	float theta = (t - (t0 + dt*idx))/((float)dt);
	//printf("hi %f %d\n", theta, idx);
	return alts[idx] + theta * (alts[idx+1] - alts[idx]);
}

template<class Float>
PressureTable<Float>::PressureTable(const char *fname) {
	FILE *f = fopen(fname, "rb");
	assert(f != 0);

	Assert(fread(&dt, 4, 1, f) == 1);
	Assert(fread(&t0, 4, 1, f) == 1);
	Assert(fread(&n, 4, 1, f) == 1);

	alts = (float*)malloc(sizeof(float)*n);
	Assert(fread(alts, sizeof(float), n, f) == n);

	printf("Loaded %d pressures from %s, dt: %d s, t0: %d.\n", n, fname, dt, t0);
}

template<class Float>
Float PressureTable<Float>::get_pressure(int t) {
	unsigned int idx = (t-t0)/dt;
	assert(idx >= 0 && idx < n);
	return alts[idx];
}

template<class Float>
wind_vector<Float> NearestNeighborWind<Float>::get_wind(int t, Float lat, Float lon, Float pres) {
	while (files[cur_file+1].time < t) {
		cur_file++;
	}
	//printf("%d cur file %d %d\n", t, cur_file, ((int)(t-files[cur_file].time)));

	/* Get pressure level. */
	int i;
	for (i=0; i<NUM_LEVELS; i++) {
		if (pres <= LEVELS[i]) break;
	}
	assert(i > 0 && i < NUM_LEVELS);
	Float theta_pr = (pres-LEVELS[i-1])/(LEVELS[i] - LEVELS[i-1]);

	#define INTERP_ALT(dst, src, idx) \
		Float dst = src->data[i-1][idx]/100. + theta_pr*(src->data[i][idx]/100. - src->data[i-1][idx]/100.);
	#define LAT(x) (LAT_MIN + LAT_D * (x))
	#define LON(x) (LON_MIN + LON_D * (x))

	point pt = get_base_neighbor(VAL(lat), VAL(lon));
	Float theta_lat = (lat - (LAT_MIN + LAT_D * pt.lat))/LAT_D;
	Float theta_lon = (lon - (LON_MIN + LON_D * pt.lon))/LON_D;
	float theta_t = t-files[cur_file].time;
	theta_t /= (files[cur_file+1].time-files[cur_file].time);
	assert(theta_t >= 0 && theta_t <= 1);

	data_file *f = files + cur_file;
	Float us[2];
	Float vs[2];
	for (int j=0; j<2; j++) {
		wind_data *p11 = get_data_at_point(f+j, {pt.lat, pt.lon});
		wind_data *p12 = get_data_at_point(f+j, {pt.lat, pt.lon+1});
		wind_data *p21 = get_data_at_point(f+j, {pt.lat+1, pt.lon});
		wind_data *p22 = get_data_at_point(f+j, {pt.lat+1, pt.lon+1});

		INTERP_ALT(u11, p11, 0);
		INTERP_ALT(v11, p11, 1);

		INTERP_ALT(u12, p12, 0);
		INTERP_ALT(v12, p12, 1);

		INTERP_ALT(u21, p21, 0);
		INTERP_ALT(v21, p21, 1);

		INTERP_ALT(u22, p22, 0);
		INTERP_ALT(v22, p22, 1);

		Float ulat1 = u11 + theta_lat * (u21 - u11);
		Float ulat2 = u12 + theta_lat * (u22 - u12);
		us[j] = ulat1 + theta_lon * (ulat2 - ulat1);

		Float vlat1 = v11 + theta_lat * (v21 - v11);
		Float vlat2 = v12 + theta_lat * (v22 - v12);
		vs[j] = vlat1 + theta_lon * (vlat2 - vlat1);
	}

	/*printf("u: %.2f %.2f\n   %.2f %.2f\n", u11, u12, u21, u22);
	printf("a: %.1f o: %.1f %.1f\n   %.1f\n", LAT(pt.lat), LON(pt.lon), LON(pt.lon+1), LAT(pt.lat+1));
	*/
		//printf(" --> %.2f @ (%f, %f)\n", u, lat, lon);
	
	//printf("theta lat %f %f [%f|%f] near %f\n", theta_lat, lat, LAT_MIN + LAT_D * pt.lat, LAT_MIN + LAT_D * (pt.lat+1), LAT_MIN + LAT_D * npt.lat);
	//printf("theta lon %f %f [%f|%f] near %f\n", theta_lon, lon, LON_MIN + LON_D * pt.lon, LON_MIN + LON_D * (pt.lon+1), LON_MIN + LON_D * npt.lon);
	
	//wind_data *data = get_data_at_point(files + cur_file, pt);

	//printf("got pt %d %d file %d\n", pt.lat, pt.lon, cur_file);
	wind_vector<Float> w;
	//w.u = data->data[i][0]/100.;
	//w.v = data->data[i][1]/100.;
	//w.u = data->data[i-1][0]/100. + theta_pr*(data->data[i][0]/100. - data->data[i-1][0]/100.);
	//w.v = data->data[i-1][1]/100. + theta_pr*(data->data[i][1]/100. - data->data[i-1][1]/100.);
	w.u = us[0] + theta_t * (us[1] - us[0]);
	w.v = vs[0] + theta_t * (vs[1] - vs[0]); 
	
	//printf("sides %f %f, %f --> %f\n", us[0], us[1], theta_t, w.u);
	//printf("hmmmm %f [%f|%f]:%f [%f|%f] -> %f\n", pres, LEVELS[i-1], LEVELS[i], theta_pr, data->data[i-1][0]/100., data->data[i][0]/100., w.u);
	return w;
}

template<class Float>
Simulation<Float>::Simulation(PressureSource<Float>& s, WindSource<Float>& w, int i)
		: pressure(s), winds(w) {
	char path[PATH_MAX];
	snprintf(path, PATH_MAX, "../ignored/output.%d.bin", i);
	file = fopen(path, "wb");
	assert(file != 0);
}

template<class Float>
vec2<Float> Simulation<Float>::run(int t, Float lat, Float lon) {
	int Tmax = tmax + t;
	const float idlat = dt / (2 * M_PI * 6371008 / 360.);
	//printf("Starting from (%f, %f)\n", VAL(lat), VAL(lon));
	clock_t t0 = clock();
	while (t < Tmax) {
		float actual_lat = VAL(lat);
		float actual_lon = VAL(lon);
		fwrite(&actual_lat, sizeof(float), 1, file);
		fwrite(&actual_lon, sizeof(float), 1, file);
		Float p = pressure.get_pressure(t);


		float actual_p = VAL(p);
		fwrite(&actual_p, sizeof(float), 1, file);

		// jank hack for simmed lasagna here
		
		//printf("lat %f lon %f pres %f\n", lat, lon, p);
		wind_vector<Float> w = winds.get_wind(t, lat, lon, p);
		lat += w.v * idlat;
		lon += w.u * idlat / cos(lat * M_PI / 180.);
		t += dt;
	}
	float dt = (clock() - t0)/((double)CLOCKS_PER_SEC)*1000;
	(void)dt;
	//printf("Ended up in (%f, %f) after %.2f ms\n", VAL(lat), VAL(lon), dt);
	vec2<Float> ret = {lat, lon};
	fclose(file);
	return ret;
}

float p2alt(float p){
  return (1.0-(pow(((float)p/101350.0),0.190284)))*145366.45*0.3048;
}

template class PressureTable<float>;
template class WaypointController<float>;
template class NearestNeighborWind<float>;
template class Simulation<float>;


#define INIT(type) \
		template class PressureTable<type>; \
		template class WaypointController<type>; \
		template class NearestNeighborWind<type>; \
		template class Simulation<type>;

INIT(adouble)
INIT(float)

