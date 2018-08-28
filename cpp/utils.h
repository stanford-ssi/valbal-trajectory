#ifndef UTILS_H
#define UTILS_h

#ifdef DEBUG_PRINT
	#define debugf(format, etc...) printf(format, ##etc);	
#else
	#define debugf(format, etc...)
#endif

#endif
