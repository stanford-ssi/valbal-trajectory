#include <stdio.h>

#include "data.h"
#include "sim.h"

using namespace std;

int main() {
	printf("ValBal trajectory optimization.\n");
	printf("This program is highly cunning and, on the face of it, not entirely wrong.\n");
	load_data();
	printf("%p\n", files);
	wind_data *sample = get_data_at_point(files+0, {42, 120});
	printf("(u,v) %hd %hd\n", sample->data[4][0], sample->data[4][1]);

	point base = get_base_neighbor(69.5, 60.6);
	point near = get_nearest_neighbor(69.5, 60.9);
	printf("(%d,%d) (%d, %d)\n", base.lat, base.lon, near.lat, near.lon);

	AltitudeTable<float> alts("../proc/alts.bin");
	Simulation<float> sim(alts);
	printf("getting altitude %f\n", sim.altitude.get_altitude(3));
}
