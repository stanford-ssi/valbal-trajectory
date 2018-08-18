#ifndef SIM_H
#define SIM_H

#include <stdint.h>

template <class Float>
class AltitudeSource {
public:
	virtual Float get_altitude(uint64_t) = 0;
};

template <class Float>
class WindSource {
public:

};

template <class Float>
class Simulation {
public:
	Simulation(AltitudeSource<Float>& s) : altitude(s) {};
	AltitudeSource<Float>& altitude;
	WindSource<Float> winds;

	int cur_file;

};

template <class Float>
class AltitudeTable : public AltitudeSource<Float> {
public:
	AltitudeTable(const char *);
	Float get_altitude(uint64_t);

private:
	uint32_t dt;
	uint32_t t0;
	uint32_t n;
	float *alts;
};

template <class Float>
class NearestNeighborWind : public WindSource<Float> {
public:

};

#endif
