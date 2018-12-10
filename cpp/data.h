#ifndef DATA_H
#define DATA_H

#include <stdint.h>

//#include "../ignored/proc/gfs_anl_0deg5/gfs_anl_0deg5.h"
//#include "../ignored/proc/gfs_anl_0deg5/gfs_anl_0deg5.h"
//#include "../ignored/proc/euro_anl/euro_anl.h"
//#include "../ignored/proc/euro_fc/euro_fc.h"
#include "config.h"

typedef short wind_t;
const float LON_MIN = 0.000000;
const float LON_MAX = 359.500000;
const float LON_D = 0.500000;
const int NUM_LONS = 720;
const float LAT_MIN = 90.000000;
const float LAT_MAX = -90.000000;
const float LAT_D = -0.500000;
const int NUM_LATS = 361;
const float LEVELS[] = {2000,3000,5000,7000,10000,15000,20000,25000,30000,35000,40000,45000,50000,55000,60000,65000,70000,75000,80000,85000,90000,92500,95000,97500,100000}; 
const int NUM_LEVELS = sizeof(LEVELS)/sizeof(LEVELS[0]); 
const int NUM_VARIABLES = 2; 

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
struct sim_state {
	Float lat;
	Float lon;
	Float p;
	Float bal;
	Float bal_rate;
	int t;
};

void load_data(const char *dname, uint64_t start, uint64_t end);
wind_t *get_data_at_point(data_file*, point);
point get_base_neighbor(float, float);
point get_nearest_neighbor(float, float);

extern data_file *files;
extern int num_files;

#endif
