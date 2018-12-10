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
#include <iostream>
#include <fstream>
#include <jsoncpp/json/json.h>
#include <sstream>
#include <algorithm>
#include <vector>
#include "data.h"


void DataHandler::load_data(const char *dname, uint64_t start, uint64_t end) {
	assert(!loaded && "bruh don't load data more than once with the same object");
	loaded = true;
	#ifdef SHOULD_PRELOAD
		struct rlimit rl;
		rl.rlim_cur = rl.rlim_max = 1024*1024*512;
		assert(setrlimit(RLIMIT_MEMLOCK, &rl) == 0 && "Need root to increase memory limit");
	#endif
	DIR *dir = opendir(dname);
	assert(dir != 0);

	char confp[PATH_MAX];
	snprintf(confp, PATH_MAX, "%s/config.json", dname);
	std::ifstream configf(confp);
	assert(!configf.fail());
	Json::Reader reader;
	Json::Value obj;
	reader.parse(configf,obj);
	std::istringstream converter(obj["hash"].asString());
	converter >> std::hex >> conf_hash; 

	LON_MIN =   obj["LON_MIN"].asFloat();		
	LON_MAX =   obj["LON_MAX"].asFloat();		
	LON_D =   	obj["LON_D"].asFloat();	
	NUM_LONS =  obj["NUM_LONS"].asInt(); 		
	LAT_MIN =   obj["LAT_MIN"].asFloat();		
	LAT_MAX =   obj["LAT_MAX"].asFloat();		
	LAT_D =   	obj["LAT_D"].asFloat();	
	NUM_LATS =  obj["NUM_LATS"].asInt(); 		
	NUM_LEVELS = obj["NUM_LEVELS"].asInt();
	LEVELS = new float[NUM_LEVELS];
	const char* l1 = obj["LEVELS"].asCString(); l1 = l1 +1;
	for(int i = 0;i<NUM_LEVELS-1;i++){
		LEVELS[i] = atof(l1);
		const char* l2 = strchr(l1,',');
		l1 = l2+1;
	}
	LEVELS[NUM_LEVELS-1] = atof(l1);

	NUM_VARIABLES = obj["NUM_VARIABLES"].asInt();

	std::vector<uint64_t> timestamps;
	struct dirent *entry;	
	while ((entry = readdir(dir)) != 0) {
		if (entry->d_name[0] < '0' || entry->d_name[0] > '9') {
			continue;
		}
		uint64_t ftime = atoll(entry->d_name);
		if (start <= ftime && ftime <= end ){
			timestamps.push_back(atoll(entry->d_name));
		}
	}
	std::sort(timestamps.begin(), timestamps.end());

	files = (data_file*)malloc(timestamps.size() * sizeof(data_file));

	for (size_t i=0; i<timestamps.size(); i++) {
		uint64_t time = timestamps[i];
		char path[PATH_MAX];
		snprintf(path, PATH_MAX, "%s/%lu.bin", dname, time);

		files[i].time = time;

		files[i].fd = open(path, O_RDONLY);
		assert(files[i].fd >= 0);

		struct stat s;
		fstat(files[i].fd, &s);
		int mmap_flags = MAP_SHARED;
		#if __linux__
			mmap_flags |= MAP_POPULATE; /* Possibly only if SHOULD_PRELOAD? */
		#endif
		files[i].data = (wind_t*)((char*)mmap(NULL, s.st_size, PROT_READ, mmap_flags, files[i].fd, 0) + 4);
		assert(files[i].data != MAP_FAILED);
		#ifdef SHOULD_PRELOAD
			assert(mlock(files[i].data, s.st_size) == 0);
		#endif
	}
	num_files = timestamps.size();
	printf("Loaded %d wind data files.\n", num_files);
}

wind_t* DataHandler::get_data_at_point(data_file *file, point p) {
	//printf("%d %d %lu %p\n", latidx, lonidx, sizeof(wind_data), file);
	if (!file->verified){
		//printf("checking file\n");
		uint32_t fhash = *(uint32_t*)((char*)file->data-4);
		fhash = (fhash & 0x000000FFU) << 24 | (fhash & 0x0000FF00U) << 8 | (fhash & 0x00FF0000U) >> 8 | (fhash & 0xFF000000U) >> 24; //god damn it
		//printf("[file hash] %x,%x\n",fhash,conf_hash);
		assert(fhash == conf_hash);
		file->verified = true;
	}
	if (p.lon >= NUM_LONS) p.lon -= NUM_LONS;
	if (p.lon < 0) p.lon += NUM_LONS;
	if (p.lat >= NUM_LATS) p.lat -= NUM_LATS;
	if (p.lat < 0) p.lat += NUM_LATS;
	return &file->data[NUM_LEVELS*NUM_VARIABLES*(NUM_LONS * p.lat + p.lon)];
}

point DataHandler::get_base_neighbor(float lat, float lon) {
	int lat0 = (int)((lat - LAT_MIN)/LAT_D);
	int lon0 = (int)((lon - LON_MIN)/LON_D);
	return {lat0, lon0};
}

point DataHandler::get_nearest_neighbor(float lat, float lon) {
	int lat0 = (int)round((lat - LAT_MIN)/LAT_D);
	int lon0 = (int)round((lon - LON_MIN)/LON_D);
	return {lat0, lon0};
}
