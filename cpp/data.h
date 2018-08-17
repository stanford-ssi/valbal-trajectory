#ifndef DATA_H
#define DATA_H

#include <stdint.h>

typedef struct {
	uint64_t time;
	short *data;
	int fd;
} data_file_t;

void load_data();

#endif
