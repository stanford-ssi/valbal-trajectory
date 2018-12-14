#include "opt.h"

void GradStep::new_step() {
	lr_set *= params.alpha;
	lr_tol *= params.alpha;
}

void GradStep::optimize(ctrl_cmd<adept::adouble>* cmd, int N) {

	for(int i = 0; i < N; i++){
		double grad_h = cmd[i].h.get_gradient();
		double grad_tol = cmd[i].tol.get_gradient();

		double val = VAL(cmd[i].h) + lr_set * grad_h;
		val = min(params.set_max, max(params.set_min, val));
		cmd[i].h.set_value(val);

		double tval = VAL(cmd[i].tol) + lr_tol * grad_tol;
		tval = min(params.tol_max, max(params.tol_min, tval));
		cmd[i].tol.set_value(tval);
	}
	//printf("[cmds] setpoint:%f, tol:%f\n",val,tval);
}
