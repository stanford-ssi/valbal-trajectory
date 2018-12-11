#ifndef TRAJ_UTILS_H
#define TRAJ_UTILS_H

#include <time.h>
#include <cmath>
#include <thread>
#include <utility>
#include <mutex>
#include <queue>
#include <functional>
#include <condition_variable>
#include "trajtypes.h"

using namespace std;

inline float VAL(float f) { return f; }
inline float VAL(adouble f) { return f.value(); }
inline adouble cosf(adouble f) { return cos(f); }
inline float fastcos(float f) { return __builtin_cos(f); }
inline adouble fastcos(adouble f) { return cos(f); }

#ifdef DEBUG_PRINT
	#define debugf(format, etc...) printf(format, ##etc);	
#else
	#define debugf(format, etc...)
#endif

/* Allows for composite variable names in macros. */
#define _JOIN(X,Y) X##Y
#define JOIN(X,Y) _JOIN(X,Y)

/* Timing utilities. TIMEIT(name, block) can wrap around a block and give timing information for
 * it. walltime and cputime return wall and CPU time, respectively, in microseconds from an
 * arbitrary and irrelevant origin. */
inline long walltime() {
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
        return (ts.tv_sec * 1000000 + ts.tv_nsec / 1000);
    } else return 0;
}

inline long cputime() {
	return clock()/((double)CLOCKS_PER_SEC)*1000000;
}

#define TIMEIT(s, blk) \
	long JOIN(cputimer,__LINE__) = cputime(); \
	long JOIN(walltimer,__LINE__) = walltime(); \
	blk; \
	printf("[%s] CPU %.2f ms, wall %.2f ms.\n", s, \
			(cputime() - JOIN(cputimer,__LINE__))/1000., \
			(walltime() - JOIN(walltimer,__LINE__))/1000.);

/* This macro is like assert, but when program logic is included inside the statement. This ensures
 * that compiling with -DNDEBUG won't make it go away. Useful for making sure system calls have the
 * proper return value. */
#define ensure(x) { \
	if (!(x)) { \
		printf("fatal error on file %s, line %d\n", __FILE__, __LINE__); \
		exit(1); \
	} \
}

/* Returns the full path to the ignored data folder. Implemented as a macro to avoid having to
 * allocate a string in the heap and keep track of it at runtime. */
#define get_data_path(x) ("../ignored/" x)

template<class ReturnType>
class Scheduler {
public:
	Scheduler(float num, int num_outputs=1024) {
		if (num < 0) {
			unsigned int logical_threads = thread::hardware_concurrency() >> 1;
			num_threads = round(logical_threads * -num);
		} else {
			num_threads = num;
		}
		results = (ReturnType*)malloc(sizeof(ReturnType)*num_outputs);
		for (int i=0; i<num_threads; i++) {
			printf("starting %d %d\n", i, num_threads);
			pool.push_back(thread([this]() {
				while (true) {
					pair<int, function<ReturnType()>> fn;
					{
						lock_guard<mutex> lg(m);
						cv.wait(m, [this] {
							return done || (tasks.size() != 0);
						});
						if (done && tasks.size() == 0) break;
						fn = tasks.front();
						tasks.pop();
					}
					results[fn.first] = fn.second();
				}
			}));
		}
	};

	void add(function<ReturnType()> fn) {
		{
			lock_guard<mutex> lg(m);
			tasks.push(make_pair(current++, fn));
		}
		cv.notify_one();
	};

	void finish() {
		done = true;
		cv.notify_all();
		for (thread& t : pool) {
			t.join();
		}
	}

	int num_threads;
	int current = 0;
	bool done = false;

	ReturnType *results;

private:
	condition_variable_any cv;
	mutex m;
	vector<thread> pool;
	queue<pair<int, function<ReturnType()>>> tasks;
};

#endif
