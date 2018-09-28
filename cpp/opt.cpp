#include "opt.h"

template<class Float>
void FinalLongitude<Float>::update(Float& lat, Float& lon, Float& p) {
	this->lon = lon;
}

template<class Float>
Float FinalLongitude<Float>::getObjective() {
	return lon;
}


template<class Float>
void MinDistanceToPoint<Float>::update(Float& lat, Float& lon, Float& p) {
	Float dist = sqrt(pow(loc[0]-lat,2.) + pow(loc[1]-lon,2.));
	min_dist = min(min_dist,dist);
	//printf("%f,%f : %f,%f\n",VAL(lat),VAL(lon),VAL(loc[0]),VAL(loc[1]));
}

template<class Float>
Float MinDistanceToPoint<Float>::getObjective() {
	return min_dist;
}


#define INIT_OPT(type) \
		template class NoOp<type>;\
		template class FinalLongitude<type>;\
		template class MinDistanceToPoint<type>;

INIT_OPT(adept::adouble)
INIT_OPT(float)