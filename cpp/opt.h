#ifndef OPT_H
#define OPT_H

#include <adept.h>
#include "trajtypes.h"
using adept::adouble;

#include "utils.h"


class StepRule {
public:
	virtual void optimize(ctrl_cmd<adept::adouble>&) = 0;
	virtual void reset() = 0;
};


class GradStep : public StepRule{
public:
	typedef struct {
		double lr_set  = 200000;
		double lr_tol  = 20000;
		double alpha   = 0.999;
		double tol_min = 200;
		double tol_max = 4000;
		double set_min = 10000;
		double set_max = 17000;
	} HyperParams;
	GradStep(){reset();};
	GradStep(HyperParams p) : params(p) {reset();};
	void reset(){lr_set = params.lr_set ; lr_tol = params.lr_tol;};
	void optimize(ctrl_cmd<adept::adouble>&);
	HyperParams params;

	double lr_set;
	double lr_tol;
};

#endif