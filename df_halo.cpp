#include "df_halo.h"
#include "math_core.h"
#include "math_ode.h"
#include <cmath>
#include <stdexcept>

#define DE

namespace df{

namespace {  // internal ns

/// helper class used in the root-finder to determine the auxiliary
/// coefficient beta for a cored halo
/// We ned one for DoublePowerLaw and one for ModifiedDoublePowerLaw
class BetaFinder: public math::IFunctionNoDeriv{
	const DoublePowerLawParam& par;

    // return the difference between the non-modified and modified DF as a function of beta and
    // the appropriately scaled action variable (t -> hJ), weighted by d(hJ)/dt for the integration in t
	double deltaf(const double t, const double beta) const
	{
	// integration is performed in a scaled variable t, ranging from 0 to 1,
	// which is remapped to hJ ranging from 0 to infinity as follows:
		double hJ    = par.Jcore * t*t*(3-2*t) / pow_2(1-t) / (1+2*t); 
		double dhJdt = par.Jcore * 6*t / pow_3(1-t) / pow_2(1+2*t);
		return hJ * hJ * dhJdt *
				math::pow(1 + math::pow(par.J0 / hJ, par.steepness),  par.slopeIn  / par.steepness) *
				math::pow(1 + math::pow(hJ / par.J0, par.steepness), -par.slopeOut / par.steepness) *
				(math::pow(1 + par.Jcore/hJ * (par.Jcore/hJ - beta), -0.5*par.slopeIn) - 1);
	}

	public:
		BetaFinder(const DoublePowerLawParam& _par) : par(_par) {}

		virtual double value(const double beta) const
		{
			double result = 0;
	// use a fixed-order GL quadrature to compute the integrated difference in normalization between
	// unmodified and core-modified DF, which is sought to be zero by choosing an appropriate beta
			static const int GLORDER = 20;  // should be even, to avoid singularity in the integrand at t=0.5
			for(int i=0; i<GLORDER; i++)
				result += math::GLWEIGHTS[GLORDER][i] * deltaf(math::GLPOINTS[GLORDER][i], beta);
			return result;
		}
};

// helper function to compute the auxiliary coefficient beta in the case of a central core
double computeBeta(const DoublePowerLawParam &par)
{
	if(par.Jcore<=0)
		return 0;
	return math::findRoot(BetaFinder(par), 0.0, 2.0, /*root-finder tolerance*/ SQRT_DBL_EPSILON);
}

/// 2nd helper class used in the root-finder to determine the auxiliary coefficient beta for a cored halo
class ModifiedBetaFinder: public math::IFunctionNoDeriv{
	const ModifiedDoublePowerLawParam& par;

    // return the difference between the non-modified and modified DF as a function of beta and
    // the appropriately scaled action variable (t -> hJ), weighted by d(hJ)/dt for the integration in t
	double deltaf(const double t, const double beta) const
	{
	// integration is performed in a scaled variable t, ranging from 0 to 1,
	// which is remapped to hJ ranging from 0 to infinity as follows:
		double hJ    = par.Jcore * t*t*(3-2*t) / pow_2(1-t) / (1+2*t); 
		double dhJdt = par.Jcore * 6*t / pow_3(1-t) / pow_2(1+2*t);
		return hJ * hJ * dhJdt *
				math::pow(1 + par.J0 / hJ,  par.slopeIn) *
				math::pow(1 + hJ / par.J0, -par.slopeOut) *
				(math::pow(1 + par.Jcore/hJ * (par.Jcore/hJ - beta), -0.5*par.slopeIn) - 1);
	}

	public:
		ModifiedBetaFinder(const ModifiedDoublePowerLawParam& _par) : par(_par) {}

		virtual double value(const double beta) const
		{
			double result = 0;
	// use a fixed-order GL quadrature to compute the integrated difference in normalization between
	// unmodified and core-modified DF, which is sought to be zero by choosing an appropriate beta
			static const int GLORDER = 20;  // should be even, to avoid singularity in the integrand at t=0.5
			for(int i=0; i<GLORDER; i++)
				result += math::GLWEIGHTS[GLORDER][i] * deltaf(math::GLPOINTS[GLORDER][i], beta);
			return result;
		}
};

// helper function to compute the auxiliary coefficient beta in the case of a central core
double computeBeta(const ModifiedDoublePowerLawParam &par)
{
	if(par.Jcore<=0)
		return 0;
	return math::findRoot(ModifiedBetaFinder(par), 0.0, 2.0, /*root-finder tolerance*/ SQRT_DBL_EPSILON);
}

}  // internal ns

EXP DoublePowerLaw::DoublePowerLaw(const DoublePowerLawParam &inparams) :
    par(inparams), beta(computeBeta(par))
{
    // sanity checks on parameters
	if(!(par.norm>0))
		throw std::invalid_argument("DoublePowerLaw: normalization must be positive");
	if(!(par.J0>0))
		throw std::invalid_argument("DoublePowerLaw: break action J0 must be positive");
	if(!(par.Jcore>=0 && beta>=0))
		throw std::invalid_argument("DoublePowerLaw: core action Jcore is invalid");
	if(!(par.Jcutoff>=0))
		throw std::invalid_argument("DoublePowerLaw: truncation action Jcutoff must be non-negative");
	if(!(par.slopeOut>3) && par.Jcutoff==0)
		throw std::invalid_argument(
					    "DoublePowerLaw: mass diverges at large J (outer slope must be > 3)");
	if(!(par.slopeIn<3))
		throw std::invalid_argument(
					    "DoublePowerLaw: mass diverges at J->0 (inner slope must be < 3)");
	if(!(par.steepness>0))
		throw std::invalid_argument("DoublePowerLaw: transition steepness parameter must be positive");
	if(!(par.cutoffStrength>0))
		throw std::invalid_argument("DoublePowerLaw: cutoff strength parameter must be positive");
	if(!(par.coefJrIn>0 && par.coefJzIn >0 && par.coefJrIn +par.coefJzIn <3 &&
	     par.coefJrOut>0 && par.coefJzOut>0 && par.coefJrOut+par.coefJzOut<3) )
		throw std::invalid_argument(
					    "DoublePowerLaw: invalid weights in the linear combination of actions");
	if(!(fabs(par.rotFrac)<=1))
		throw std::invalid_argument(
					    "DoublePowerLaw: amplitude of odd-Jphi component must be between -1 and 1");
	if(par.Fname!=""){
		readBrighterThan(par.Fname);
	}
}

EXP double DoublePowerLaw::value(const actions::Actions &J) const
{
    // linear combination of actions in the inner part of the model (for J<J0)
	double hJ  = par.coefJrIn * J.Jr + par.coefJzIn * J.Jz +
		     (3-par.coefJrIn -par.coefJzIn) * fabs(J.Jphi);
    // linear combination of actions in the outer part of the model (for J>J0)
	double gJ  = par.coefJrOut* J.Jr + par.coefJzOut* J.Jz +
		     (3-par.coefJrOut-par.coefJzOut)* fabs(J.Jphi);
	double val = par.norm / pow_3(2*M_PI * par.J0) *
		     math::pow(1 + math::pow(par.J0 / hJ, par.steepness),  par.slopeIn  / par.steepness) *
		     math::pow(1 + math::pow(gJ / par.J0, par.steepness), -par.slopeOut / par.steepness);
	if(par.rotFrac!=0)  // add the odd part
		val *= 1 + par.rotFrac * tanh(J.Jphi / par.Jphi0);
	if(par.Jcutoff>0)   // exponential cutoff at large J
		val *= exp(-math::pow(gJ / par.Jcutoff, par.cutoffStrength));
	if(par.Jcore>0) {   // central core of nearly-constant f(J) at small J
		if(hJ==0) return par.norm / pow_3(2*M_PI * par.J0);
		val *= math::pow(1 + par.Jcore/hJ * (par.Jcore/hJ - beta), -0.5*par.slopeIn);
	}
	return val;
}

EXP ModifiedDoublePowerLaw::ModifiedDoublePowerLaw(const ModifiedDoublePowerLawParam &inparams) :
    par(inparams), beta(computeBeta(par))
{
    // sanity checks on parameters
	if(!(par.norm>0))
		throw std::invalid_argument("ModifiedDoublePowerLaw: normalization must be positive");
	if(!(par.J0>0))
		throw std::invalid_argument("ModifiedDoublePowerLaw: break action J0 must be positive");
	if(!(par.Jcore>=0 && beta>=0))
		throw std::invalid_argument("ModifiedDoublePowerLaw: core action Jcore is invalid");
	if(!(par.Jcutoff>=0))
		throw std::invalid_argument("ModifiedDoublePowerLaw: truncation action Jcutoff must be non-negative");
	if(!(par.slopeOut>3) && par.Jcutoff==0)
		throw std::invalid_argument(
					    "ModifiedDoublePowerLaw: mass diverges at large J (outer slope must be > 3)");
	if(!(par.slopeIn<3))
		throw std::invalid_argument(
					    "ModifiedDoublePowerLaw: mass diverges at J->0 (inner slope must be < 3)");
	if(!(par.cutoffStrength>0))
		throw std::invalid_argument("ModifiedDoublePowerLaw: cutoff strength parameter must be positive");
	if(!(fabs(par.rotFrac)<=1))
		throw std::invalid_argument(
					    "ModifiedDoublePowerLaw: amplitude of odd-Jphi component must be between -1 and 1");
	if(par.Fname!=""){
		readBrighterThan(par.Fname);
	}
}

EXP double ModifiedDoublePowerLaw::value(const actions::Actions &J) const
{
	double L=J.Jz+fabs(J.Jphi);
	double c=L/(L+J.Jr), i=J.Jz/L;
	double jt=(1.5*J.Jr+L)/par.L0, zeta=jt/(0.5+jt);
	double jta=pow(jt,par.alpha), xi=jta/(1+jta);
	double rat=(1-zeta)*par.Fin+zeta*par.Fout;
	double a=.5*(rat+1), b=.5*(rat-1);
	double cL=J.Jz*(a+b*fabs(J.Jphi)/L)+fabs(J.Jphi);
	double fac=exp(par.beta*sin(0.5*M_PI*c));
	double hJ=J.Jr/fac + .5*(1+c*xi)*fac*cL;
	double gJ=hJ;
	double val = par.norm / pow_3(2*M_PI * par.J0) *
		     math::pow(1 + par.J0 / hJ,  par.slopeIn) *
		     math::pow(1 + gJ / par.J0, -par.slopeOut);
	if(par.Jcutoff>0)   // exponential cutoff at large J
		val *= exp(-math::pow(gJ / par.Jcutoff, par.cutoffStrength));
	if(par.Jcore>0) {   // central core of nearly-constant f(J) at small J
		if(hJ==0) return par.norm / pow_3(2*M_PI * par.J0);
		val *= math::pow(1 + par.Jcore/hJ * (par.Jcore/hJ - beta), -0.5*par.slopeIn);
	}
	if(par.rotFrac!=0)  // add the odd part
		val *= 1 + par.rotFrac * tanh(J.Jphi / par.Jphi0);
	return val;
}
}  // namespace df
