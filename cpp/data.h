#ifndef DATA_H
#define DATA_H

#include <stdint.h>

#include "../ignored/proc/gfs_anl_0deg5/gfs_anl_0deg5.h"
//#include "../ignored/proc/euro_anl/euro_anl.h"
//#include "../ignored/proc/euro_fc/euro_fc.h"
#include "config.h"

typedef struct __attribute__((packed)) {
	wind_t data[NUM_LEVELS][NUM_VARIABLES];
} wind_data;

typedef struct {
	int64_t time;
	wind_data *data;
	int fd;
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

void load_data(const char *dname, uint64_t start, uint64_t end);
wind_data *get_data_at_point(data_file*, point);
point get_base_neighbor(float, float);
point get_nearest_neighbor(float, float);

extern data_file *files;
extern int num_files;

#endif
