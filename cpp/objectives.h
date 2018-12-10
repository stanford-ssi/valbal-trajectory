#ifndef OBJECTIVES_H
#define OBJECTIVES_H

#include <random>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <adept.h>
using adept::adouble;
#include "utils.h"
#include "trajtypes.h"


template <class Float>
class ObjectiveFn {
public:
	virtual Float update(sim_state<Float>& state, bool save = true) = 0;
	virtual Float getObjective() = 0;
};

template <class Float>
class FinalLongitude : public ObjectiveFn<Float> {
public:
	FinalLongitude(){};
	Float update(sim_state<Float>& state, bool save = true);
	Float getObjective();
	Float lon = 0;
};

template <class Float>
class MinDistanceToPoint : public ObjectiveFn<Float> {
public:
	MinDistanceToPoint(float lat, float lon) : loc{lat,lon} {};
	Float update(sim_state<Float>& state, bool save = true);
	Float getObjective();
	float loc[2];
	Float min_dist = 1000000;
};

template <class Float>
class DirectionToTarget : public ObjectiveFn<Float> {
public:
	DirectionToTarget(float lat, float lon) : loc{lat,lon} {};
	Float update(sim_state<Float>& state, bool save = true);
	Float getObjective();
	float loc[2];
};

template <class Float>
class NoOp : public ObjectiveFn<Float> {
public:
	NoOp(){};
	Float update(sim_state<Float>& state, bool save = true){Float ret = 0; return ret;};
	Float getObjective(){Float ret = 0; return ret;};
};



#endif