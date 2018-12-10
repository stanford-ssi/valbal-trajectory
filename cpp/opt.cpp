#include "opt.h"


void GradStep::optimize(ctrl_cmd<adept::adouble>& cmd) {
	lr_set *= params.alpha;
	lr_tol *= params.alpha;

	double grad_h = cmd.h.get_gradient();
	double grad_tol = cmd.tol.get_gradient();

	double val = VAL(cmd.h) + lr_set * grad_h;
	val = min(params.set_max, max(params.set_min, val));
	cmd.h.set_value(val);

	double tval = VAL(cmd.tol) + lr_tol * grad_tol;
	tval = min(params.tol_max, max(params.tol_min, tval));
	cmd.tol.set_value(tval);
	//printf("[cmds] setpoint:%f, tol:%f\n",val,tval);
}
