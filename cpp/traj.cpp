#include <stdio.h>

#include "data.h"
#include "sim.h"

using namespace std;

int main() {
	printf("ValBal trajectory optimization.\n");
	printf("This program is highly cunning and, on the face of it, not entirely wrong.\n");

	load_data("../ignored/proc/gfs_anl_0deg5", 1500000000,1600000000);

	/*wind_data *sample = get_data_at_point(files+0, {42, 120});
	printf("(u,v) %hd %hd\n", sample->data[4][0], sample->data[4][1]);

	point base = get_base_neighbor(69.5, 60.6);
	point near = get_nearest_neighbor(69.5, 60.9);
	printf("(%d,%d) (%d, %d)\n", base.lat, base.lon, near.lat, near.lon);
	*/

	PressureTable<float> pres("../proc/inp.bin");
	NearestNeighborWind<float> winds;
	Simulation<float> sim(pres, winds, 0);
	//sim.run(pres.t0 + pres.dt*2500, 37.7633, 239.86015);
	//printf("actual t0 %d %d\n", pres.t0 + pres.dt*5000, pres.t0);
	//sim.run(pres.t0 + pres.dt*5000, 35.339, -115.0733+360);
	sim.run(pres.t0 + pres.dt*1000, 36.95854187, -121.84505463+360);
}
