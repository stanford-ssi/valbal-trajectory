#include <stdio.h>
#include <assert.h>

#include "sim.h"

template<class Float>
AltitudeTable<Float>::AltitudeTable(const char *fname) {
	FILE *f = fopen(fname, "rb");
	assert(f != 0);

	assert(fread(&dt, 4, 1, f) == 1);
	assert(fread(&t0, 4, 1, f) == 1);
	assert(fread(&n, 4, 1, f) == 1);
	printf("Loaded %d altitudes from %s, dt: %d s, t0: %d.\n", n, fname, dt, t0);
}

template<class Float>
Float AltitudeTable<Float>::get_altitude(uint64_t t) {
	return (float)t;
}

template class AltitudeTable<float>;
