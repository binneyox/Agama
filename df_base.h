/** \file    df_base.h
    \brief   Base class for action-based distribution functions
    \authors Eugene Vasiliev, Payel Das, Fred Thompson, James Binney
    \date    2015-2021
*/
#pragma once
#include "math_spline.h"
#include "utils.h"
#include "actions_base.h"
#include "units.h"
#include <fstream>
#include <iostream>
#include <vector>
#define EXP __declspec(dllexport)

/** Classes for dealing with action-base distribution functions */
namespace df{

/** Base class defining the action-based distribution function (DF) */
class EXP BaseDistributionFunction{
	private:
		//interps(M) returns fraction of pop brighter than M
		math::LinearInterpolator interps;
		
public:
    BaseDistributionFunction() {};
    virtual ~BaseDistributionFunction() {};

    /** Compute the total mass, i.e., integral of DF over the entire phase space.
        The actions range between 0 and +inf for Jr and Jz, and between -inf and +inf for Jphi,
        and there is an additional factor (2*pi)^3 from the integration over angles
        (we assume that DF does not depend on angles, but we still need to integrate
        over the entire 6d phase space).
        Derived classes may return an analytically computed value if available,
        and the default implementation performs multidimension integration numerically.
        \param[in]  reqRelError - relative tolerance;
        \param[in]  maxNumEval - maximum number of evaluations of DF during integration;
        \param[out] error (optional) if not NULL, store the estimate of the error;
        \param[out] numEval (optional) if not NULL, store the actual number of DF evaluations.
        \returns    total mass
    */
    virtual double totalMass(const double reqRelError=1e-6, const int maxNumEval=1e6,
        double* error=NULL, int* numEval=NULL) const;

    /** Value of distribution function for the given set of actions J */
    virtual double value(const actions::Actions &J) const=0;

    /** Number of components in the case of a multi-component DF */
    virtual unsigned int numValues() const { return 1; }

    /** Compute values of all components for the given actions */
    virtual void eval(const actions::Actions &J, double values[]) const {
	    *values = value(J);
    }

    /** Compute values of all components for the given actions */
    virtual void wthSF(const actions::Actions &J, double values[], double Bright, double Faint) const {
	    *values = value(J);
    }

    virtual void set_norm(double){std::cout << "set_norm not defined\n";}
    /** write parameters */
    virtual void write_params(std::ofstream& strm,const units::InternalUnits& conv) const {strm << "write_params not defined\n";}
    virtual void tab_params  (std::ofstream& strm,const units::InternalUnits& conv) const {strm << "tab_params not defined\n";}
    /** estimate of Omega_phi/Omega_r */
    virtual void epicycle_ratios(const actions::Actions&,double*) const{}

    /* The fraction of the populaton that lies between bright and faint
     * abs mags. Influences weights with which components contribute to
     * data
     */
    virtual double fraction(double bright,double faint) const;
    
    /* the LuminosityFunction, which is the derivative of fraction wrt
     * bright
     */
    virtual double LF(double mag) const;
    
    /* Multiply previously determined f(J) by fraction of pop
     * conributing to data 
    */
    virtual void withSF(const actions::Actions &J, double values[], double Bright, double Faint) const {
	    *values = value(J) * fraction(Bright, Faint);
    }

    virtual void withLF(const actions::Actions &J, double values[], double Mag) const {
	    *values = value(J) * LF(Mag);
    }

    /* Randomly choose from LF an abs magnitude in the given range
     */
    virtual double selectMag(const double Bright,const double Faint) const;
    
    /* Read in file giving fraction brighter than M starting at most
     * luminous M
    */
    virtual void readBrighterThan(std::string Fname);
};


/** Helper class for scaling transformations in the action space,
    mapping the entire possible range of actions into a unit cube */
class EXP BaseActionSpaceScaling {
public:
    virtual ~BaseActionSpaceScaling() {}

    /** Convert from scaled variables to the values of actions.
        \param[in]  vars  are the coordinates in the unit cube;
        \param[out] jac   if not NULL, will contain the value of jacobian of transformation;
        \return  the values of actions.
    */
    virtual actions::Actions toActions(const double vars[3], double *jac=NULL) const = 0;

    /** Convert from actions to scaled variables.
        \param[in]  acts  are the actions;
        \param[out] vars  will contain the scaled variables in the unit cube;
    */
    virtual void toScaled(const actions::Actions &acts, double vars[3]) const = 0;
};

class EXP ActionSpaceScalingTriangLog: public BaseActionSpaceScaling {
public:
    virtual actions::Actions toActions(const double vars[3], double *jac=NULL) const;
    virtual void toScaled(const actions::Actions &acts, double vars[3]) const;
};

class EXP ActionSpaceScalingRect: public BaseActionSpaceScaling {
    double scaleJm, scaleJphi;
public:
    ActionSpaceScalingRect(double scaleJm, double scaleJphi);
    virtual actions::Actions toActions(const double vars[3], double *jac=NULL) const;
    virtual void toScaled(const actions::Actions &acts, double vars[3]) const;
};


/** Compute the entropy  \f$  S = -\int d^3 J f(J) ln(f(J))  \f$.
    \param[in]  DF is the distribution function;
    \param[in]  reqRelError - relative tolerance;
    \param[in]  maxNumEval - maximum number of evaluations of DF during integration;
    \return  the value of entropy S.
*/
double EXP totalEntropy(const BaseDistributionFunction& DF,
    const double reqRelError=1e-4, const int maxNumEval=1e6);

/** Sample the distribution function in actions.
    In other words, draw N sampling points from the action space, so that the density of points 
    in the neighborhood of any point is proportional to the value of DF at this point 
    (point = triplet of actions).
    \param[in]  DF  is the distribution function;
    \param[in]  numSamples  is the required number of sampling points;
    \param[out] totalMass (optional) if not NULL, will store the Monte Carlo estimate 
    of the integral of the distribution function (i.e., the same quantity as computed by 
    BaseDistributionFunction::totalMass(), but calculated with a different method).
    \param[out] totalMassErr (optional) if not NULL, will store the error estimate of the integral.
    \returns    the array of actions sampled from the DF.
*/
EXP std::vector<actions::Actions> sampleActions(const BaseDistributionFunction& DF,
    const std::size_t numSamples, double* totalMass=NULL, double* totalMassErr=NULL);

}  // namespace df
