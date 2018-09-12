#ifndef CONFIG_H
#define CONFIG_H

/* Preload the wind data into physical RAM. It's otherwise memory mapped so it has a much smaller
 * memory footprint. */
//#define SHOULD_PRELOAD

/* Stores altitude, instead of pressure, to the summary file. It's about 16% slower. */
//#define STORE_PRESSURE

/* Sprinkles random debug statements throughout the code. Significant slowdown. */
//#define DEBUG_PRINT

/* Interpolates in time. May violate causality. */
#define INTERPOLATE_IN_TIME

#endif
