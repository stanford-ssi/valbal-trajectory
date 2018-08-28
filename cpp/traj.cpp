#include <stdio.h>

#include <adept.h>
using adept::adouble;

#include "data.h"
#include "sim.h"

using namespace std;

int main() {
	printf("ValBal trajectory optimization.\n");
	printf("This program is highly cunning and, on the face of it, not entirely wrong.\n");

	load_data("../ignored/proc/gfs_anl_0deg5", 1500000000,1600000000);
	//load_data("../proc", 1500000000,1600000000);

	/*wind_data *sample = get_data_at_point(files+0, {42, 120});
	printf("(u,v) %hd %hd\n", sample->data[4][0], sample->data[4][1]);

	point base = get_base_neighbor(69.5, 60.6);
	point near = get_nearest_neighbor(69.5, 60.9);
	printf("(%d,%d) (%d, %d)\n", base.lat, base.lon, near.lat, near.lon);
	*/

	int t0 = 1512871200;
	int dt = 3600;
	const int N_W = 101;
	const float LR = 10;
	double waypoints_val[N_W]; for (int i=0; i<N_W; i++) waypoints_val[i] = 15000;

	for (int it=0; it<100; it++) {
		clock_t timer0 = clock();
		adept::Stack stack;
		adouble waypoints[N_W];
		for (int i=0; i<N_W; i++) {
			waypoints[i] = min(23000., max(10000., waypoints_val[i]));
		}
		stack.new_recording();
	
		WaypointController<adouble> pres(t0, dt, waypoints);
		NearestNeighborWind<adouble> winds;
		Simulation<adouble> sim(pres, winds, it);
		float lon0 = -121.84505463+360;
		vec2<adouble> end = sim.run(pres.t0, 36.95854187, lon0);
		adouble cost = (end.b-lon0)*111195;
		cost.set_gradient(1.0);
		stack.compute_adjoint();
		for (int i=0; i<N_W; i++) {
			waypoints_val[i] += LR * waypoints[i].get_gradient();
			//printf("Update at %d: %f\n", i, LR * waypoints[i].get_gradient());
		}
		float dt = (clock() - timer0)/((double)CLOCKS_PER_SEC)*1000;
		printf("Took %.2f ms, got %f\n", dt, VAL(cost)/1e6);
	}
}
