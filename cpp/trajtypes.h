#ifndef TRAJTYPES_H
#define TRAJTYPES_H

typedef short wind_t;

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

template<class Float>
struct ctrl_cmd {
	Float h;
	Float tol;
};

#endif