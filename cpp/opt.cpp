

#include "opt.h"

template<class Float>
Float FinalLongitude<Float>::update(sim_state<Float>& state, bool save) {
	this->lon = state.lon;
	return lon;
}

template<class Float>
Float FinalLongitude<Float>::getObjective() {
	return lon;
}


template<class Float>
Float MinDistanceToPoint<Float>::update(sim_state<Float>& state, bool save) {
	Float dist = sqrt(pow(loc[0]-state.lat,2.) + pow(loc[1]-state.lon,2.));
	if(save) min_dist = min(min_dist,dist);
	return dist;
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