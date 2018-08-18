#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <errno.h>
#include <dirent.h>
#include <limits.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <sys/resource.h>

#include <algorithm>
#include <vector>

#include "config.h"
#include "data.h"

data_file *files;
int num_files;

void load_data() {
	#ifdef SHOULD_PRELOAD
		struct rlimit rl;
		rl.rlim_cur = rl.rlim_max = 1024*1024*512;
		assert(setrlimit(RLIMIT_MEMLOCK, &rl) == 0 && "Need root to increase memory limit");
	#endif
	DIR *dir = opendir("../proc");
	assert(dir != 0);

	std::vector<uint64_t> timestamps;
	struct dirent *entry;	
	while ((entry = readdir(dir)) != 0) {
		if (entry->d_name[0] < '0' || entry->d_name[0] > '9') {
			continue;
		}
		timestamps.push_back(atoll(entry->d_name));
	}
	std::sort(timestamps.begin(), timestamps.end());

	files = (data_file*)malloc(timestamps.size() * sizeof(data_file));

	for (size_t i=0; i<timestamps.size(); i++) {
		uint64_t time = timestamps[i];
		char path[PATH_MAX];
		snprintf(path, PATH_MAX, "../proc/%lu.bin", time);

		files[i].time = time;

		files[i].fd = open(path, O_RDONLY);
		assert(files[i].fd >= 0);

		struct stat s;
		fstat(files[i].fd, &s);
		int mmap_flags = MAP_SHARED;
		mmap_flags |= MAP_POPULATE; /* Possibly only if SHOULD_PRELOAD? */
		files[i].data = (wind_data*)mmap(NULL, s.st_size, PROT_READ, mmap_flags, files[i].fd, 0);
		assert(files[i].data != MAP_FAILED);
		#ifdef SHOULD_PRELOAD
			assert(mlock(files[i].data, s.st_size) == 0);
		#endif
	}
	num_files = timestamps.size();
	printf("Loaded %d wind data files.\n", num_files);
}

wind_data *get_data_at_point(data_file *file, point p) {
	//printf("%d %d %lu %p\n", latidx, lonidx, sizeof(wind_data), file);
	return &file->data[NUM_LONS * p.lat + p.lon];
}

point get_base_neighbor(float lat, float lon) {
	int lat0 = (int)((lat - LAT_MIN)/LAT_D);
	int lon0 = (int)((lon - LON_MIN)/LON_D);
	return {lat0, lon0};
}

point get_nearest_neighbor(float lat, float lon) {
	int lat0 = (int)round((lat - LAT_MIN)/LAT_D);
	int lon0 = (int)round((lon - LON_MIN)/LON_D);
	return {lat0, lon0};
}
