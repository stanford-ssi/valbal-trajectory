#ifndef DATA_H
#define DATA_H

#include <stdint.h>

//#include "../ignored/proc/gfs_anl_0deg5/gfs_anl_0deg5.h"
//#include "../ignored/proc/gfs_anl_0deg5/gfs_anl_0deg5.h"
//#include "../ignored/proc/euro_anl/euro_anl.h"
//#include "../ignored/proc/euro_fc/euro_fc.h"
#include "config.h"
#include "trajtypes.h"

void getRecentDir(char* buf, const char* dname, int64_t time);

class DataHandler {
public:
	void load_data(const char *dname, uint64_t start, uint64_t end);
	wind_t *get_data_at_point(data_file*, point);
	point get_base_neighbor(float, float);
	point get_nearest_neighbor(float, float);
	~DataHandler();
	data_file *files;
	int num_files;
	uint32_t conf_hash;
	bool loaded = false;

	float LON_MIN = 0.000000;
	float LON_MAX = 359.500000;
	float LON_D = 0.500000;
	int NUM_LONS = 720;
	float LAT_MIN = 90.000000;
	float LAT_MAX = -90.000000;
	float LAT_D = -0.500000;
	int NUM_LATS = 361;
	float* LEVELS;
	int NUM_LEVELS;
	int NUM_VARIABLES = 2;
};
#endif
