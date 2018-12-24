#ifndef TRAJTYPES_H
#define TRAJTYPES_H

#include "utils.h"

typedef short wind_t;

typedef struct {	
	int64_t time;
	wind_t* data;
	int fd;
	bool verified = false;
} data_file;

typedef struct {
	int lat;
	int lon;
} point;

template<class Float>
struct wind_vector {
	Float u;
	Float v;
};

template<class Float>
struct vec2 {
	Float a;
	Float b;
};

template<class Float>
struct ctrl_cmd {
	Float h;
	Float tol;
	template <typename NewType> 
	ctrl_cmd<NewType> cast() const{
		ctrl_cmd<NewType> ret;
		ret.h = VAL(h);
		ret.tol = VAL(tol);
        return ret;
	}
};

template<class Float>
struct sim_state {
	Float lat;
	Float lon;
	Float p;
	Float bal;
	Float bal_rate;
	ctrl_cmd<Float> cmd;
	int t;
	template <typename NewType> 
	sim_state<NewType> cast() const{
        sim_state<NewType> ret;
        ret.lat = VAL(lat);
        ret.lon = VAL(lon);
        ret.p = VAL(p);
        ret.bal = VAL(bal);
        ret.bal_rate = VAL(bal_rate);
        ret.cmd = cmd.template cast<NewType>();
        ret.t = t;
        return ret;
    } 
};

#endif