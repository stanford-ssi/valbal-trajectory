#ifndef CONFIG_H
#define CONFIG_H

/* Preload the wind data into physical RAM. It's otherwise memory mapped so it has a much smaller
 * memory footprint. */
//#define SHOULD_PRELOAD

/* Type used to store wind data on disk, in cm/s. */
typedef short wind_t;

/* What the grid is assumed to look like. */
const float LAT_MIN = 90;
const float LAT_MAX = -90;
const float LAT_D = -0.5;
const int NUM_LATS = 361;

const float LON_MIN = 0;
const float LON_MAX = 359.5;
const float LON_D = 0.5;
const int NUM_LONS = 720;

const float LEVELS[] = {5000, 7000, 10000, 15000, 20000, 25000, 30000, 35000};

const int NUM_LEVELS = sizeof(LEVELS)/sizeof(LEVELS[0]);
const int NUM_VARIABLES = 2;

#endif
