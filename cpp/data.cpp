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

data_file_t *files;
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
	int i = 0;
	while ((entry = readdir(dir)) != 0) {
		if (entry->d_name[0] == '.') {
			continue;
		}
		timestamps.push_back(atoll(entry->d_name));
	}
	std::sort(timestamps.begin(), timestamps.end());

	files = (data_file_t*)malloc(timestamps.size() * sizeof(data_file_t));

	for (uint64_t time : timestamps) {
		char path[PATH_MAX];
		snprintf(path, PATH_MAX, "../proc/%lu.bin", time);
		printf("%s\n", path);

		// Resize the list to fit. Usually a no-op
		files = (data_file_t*)realloc(files, (i+1)*sizeof(data_file_t));
		assert(files != 0);

		files[i].time = time;

		files[i].fd = open(path, O_RDONLY);
		assert(files[i].fd >= 0);

		struct stat s;
		fstat(files[i].fd, &s);
		int mmap_flags = MAP_SHARED;
		mmap_flags |= MAP_POPULATE; /* Possibly only if SHOULD_PRELOAD? */
		files[i].data = (wind_t*)mmap(NULL, s.st_size, PROT_READ, mmap_flags, files[i].fd, 0);
		assert(files[i].data != MAP_FAILED);
		#ifdef SHOULD_PRELOAD
			assert(mlock(files[i].data, s.st_size) == 0);
		#endif

		i++;
	}
	num_files = i;
}
