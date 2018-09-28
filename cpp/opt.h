#ifndef OPT_H
#define OPT_H

#include <random>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <adept.h>
using adept::adouble;
#include "utils.h"
#include "data.h"	


template <class Float>
class ObjectiveFn {
public:
	virtual void update(Float&, Float&, Float&) = 0;
	virtual Float getObjective() = 0;
};

template <class Float>
class FinalLongitude : public ObjectiveFn<Float> {
public:
	FinalLongitude(){};
	void update(Float&, Float&, Float&);
	Float getObjective();
	Float lon = 0;
};

template <class Float>
class MinDistanceToPoint : public ObjectiveFn<Float> {
public:
	MinDistanceToPoint(float lat, float lon) : loc{lat,lon} {};
	void update(Float&, Float&, Float&);
	Float getObjective();
	float loc[2];
	Float min_dist = 1000000;
};

template <class Float>
class NoOp : public ObjectiveFn<Float> {
public:
	NoOp(){};
	void update(Float&, Float&, Float&){};
	Float getObjective(){Float ret = 0; return ret;};
};



#endif