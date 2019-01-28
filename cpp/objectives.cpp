#include "objectives.h"


template<class Float>
ObjectiveFn<Float>::~ObjectiveFn<Float>(){}

template<class Float>
Float FinalLongitude<Float>::update(sim_state<Float>& state, bool save) {
	this->lon = state.lon;
	return -lon;
}

template<class Float>
Float FinalLongitude<Float>::getObjective() {
	return -lon;
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

template <class Float>
ObjectiveFn<Float>& objParse(const std::string& conf){
	Json::Reader reader;
    Json::Value vals;
    reader.parse(conf, vals);
    string type = vals["type"].asString();
    ObjectiveFn<Float>* objfn;
    if(type.compare("MinDistanceToPoint")==0){
    	assert(vals["lon"].isDouble() && vals["lat"].isDouble() && "Must have valid lat lon for MinDistanceToPoint objective");
    	objfn = new MinDistanceToPoint<Float>(vals["lat"].asDouble(),vals["lon"].asDouble());
    } else if(type.compare("FinalLongitude")==0) {
    	objfn = new FinalLongitude<Float>;
    } else { 
    	objfn = new FinalLongitude<Float>;
    }
    return *objfn;
}

#define INIT_OPT(type) \
		template class ObjectiveFn<type>;\
		template class NoOp<type>;\
		template class FinalLongitude<type>;\
		template class MinDistanceToPoint<type>;\
		template ObjectiveFn<type>& objParse(const std::string&);\

INIT_OPT(adept::adouble)
INIT_OPT(float)