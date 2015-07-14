// An interface to GSL - hides all the gsl calls behind easier to use functions
// Jason Sanders
#if 0
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_permute.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_siman.h>
#include <vector>
#endif

#include "mathutils.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_odeiv2.h>
#include <stdexcept>
#include <cassert>

namespace mathutils{

void exceptionally_awesome_gsl_error_handler (const char *reason, const char * /*file*/, int /*line*/, int gsl_errno)
{
    if( // list error codes that are non-critical and don't need to be reported
        gsl_errno == GSL_ETOL ||
        gsl_errno == GSL_EROUND ||
        gsl_errno == GSL_ESING ||
        gsl_errno == GSL_EDIVERGE )
        return;  // do nothing
    if( gsl_errno == GSL_ERANGE ||
        gsl_errno == GSL_EOVRFLW )
        throw std::range_error(std::string("GSL range error: ")+reason);
    if( gsl_errno == GSL_EDOM )
        throw std::domain_error(std::string("GSL domain error: ")+reason);
    if( gsl_errno == GSL_EINVAL )
        throw std::invalid_argument(std::string("GSL invalid argument error: ")+reason);
    throw std::runtime_error(std::string("GSL error: ")+reason);
}

// a static variable that initializes our error handler
bool error_handler_set = gsl_set_error_handler(&exceptionally_awesome_gsl_error_handler);

bool isFinite(double x) {
    return gsl_finite(x);
}

int fcmp(double x, double y, double eps) {
    return gsl_fcmp(x, y, eps);
}

double wrapAngle(double x) {
    return gsl_sf_angle_restrict_pos(x);
}

double unwrapAngle(double x, double xprev) {
    double diff=(x-xprev)/(2*M_PI);
    double nwraps=0;
    if(diff>0.5) 
        modf(diff+0.5, &nwraps);
    else if(diff<-0.5) 
        modf(diff-0.5, &nwraps);
    return x - 2*M_PI * nwraps;
}

struct RootFinderParam {
    function fnc;
    void* param;
    double x_edge, x_scaling;
    bool inf_lower, inf_upper;
};

double scaledArgumentFnc(double y, void* v_param)
{
    RootFinderParam* param = static_cast<RootFinderParam*>(v_param);
    double x = param->inf_upper ? (param->inf_lower ? 
        param->x_scaling*(1/(1-y)-1/y) :            // (-inf,inf)
        param->x_edge + y/(1-y)*param->x_scaling) : // [x_edge, inf)
        param->x_edge - param->x_scaling*(1-y)/y;   // (-inf, x_edge]
    return param->fnc(x, param->param);
}

double findRoot(function fnc, void* params, double xlower, double xupper, double reltoler)
{
    if(reltoler<=0)
        throw std::invalid_argument("findRoot: relative tolerance must be positive");
    if(xlower>=xupper)
        throw std::invalid_argument("findRoot: invalid interval (xlower>=xupper)");
    gsl_function F;
    RootFinderParam par;
    par.inf_lower = xlower==gsl_neginf();
    par.inf_upper = xupper==gsl_posinf();
    if(par.inf_upper || par.inf_lower) {  // apply internal scaling procedure
        F.function = &scaledArgumentFnc;
        F.params = &par;
        par.fnc = fnc;
        par.param = params;
        if(par.inf_upper && !par.inf_lower) {
            par.x_edge = xlower;
            par.x_scaling = fmax(xlower*2, 1.);  // quite an arbitrary choice
        } else if(par.inf_lower && !par.inf_upper) {
            par.x_edge = xupper;
            par.x_scaling = fmax(-xupper*2, 1.);
        } else
            par.x_scaling=1;
        xlower = 0;
        xupper = 1;
    } else {  // no scaling
        F.function = fnc;
        F.params = params;
    }
    gsl_root_fsolver *solv = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    try{
        gsl_root_fsolver_set(solv, &F, xlower, xupper);
    }
    catch(std::invalid_argument&) {  // endpoints do not enclose root, or function is infinite there
        return gsl_nan();
    }
    int status=0, iter=0;
    double abstoler=fabs(xlower-xupper)*reltoler;
    //double absftoler=fmax(fabs(f1), fabs(f2)) * reltoler;
    double xroot=(xlower+xupper)/2;
    do{
        iter++;
        status = gsl_root_fsolver_iterate(solv);
        if(status!=GSL_SUCCESS) break;
        xroot  = gsl_root_fsolver_root(solv);
        xlower = gsl_root_fsolver_x_lower(solv);
        xupper = gsl_root_fsolver_x_upper(solv);
        status = gsl_root_test_interval (xlower, xupper, abstoler, reltoler);
    }
    while((status == GSL_CONTINUE /*|| fabs(fnc(xroot, params))>absftoler*/) && iter < 50);
    gsl_root_fsolver_free(solv);
    if(par.inf_upper) {
        if(par.inf_lower)
            xroot = par.x_scaling*(1/(1-xroot)-1/xroot);
        else
            xroot = par.x_edge + xroot/(1-xroot)*par.x_scaling;
    } else if(par.inf_lower)
        xroot = par.x_edge - par.x_scaling*(1-xroot)/xroot;
    return xroot;
}

double findRootGuess(function fnc, void* params, double x1, double x2, 
    double xinit, bool increasing, double reltoler)
{
    double finit=fnc(xinit, params);
    if(finit==0) 
        return xinit;
    if(!gsl_finite(finit))
        throw std::invalid_argument("findRootGuess: initial function value is not finite");
    if(x1>=x2)
        throw std::invalid_argument("findRootGuess: interval is non-positive");
    if(xinit<x1 || xinit>x2)
        throw std::invalid_argument("findRootGuess: initial guess is outside interval");
    bool searchleft = (finit<0) ^ increasing;
    if(searchleft) {
        // map the interval x:(x1,xinit) onto y:(-1,0) as  y=1/(-1 + 1/(x-xinit) - 1/(x1-xinit) )
        double y=-0.5, x, fy;
        int niter=0;
        do{
            x=1/(1/y+1/(x1-xinit)+1)+xinit;
            fy=fnc(x, params);
            if(fy==0)
                return x;
            niter++;
            if(fy*finit>0)  // go further right
                y = (y-1)/2;
        } while(fy*finit>0 && niter<60);
        if(fy*finit>0)
            return gsl_nan();
        return findRoot(fnc, params, x, xinit, reltoler);
    } else {  // search right
        // map the interval x:(xinit,x2) onto y:(0,1) as  y=1/(1 + 1/(x-xinit) - 1/(x2-xinit) )
        double y=0.5, x, fy;
        int niter=0;
        do{
            x=1/(1/y+1/(x2-xinit)-1)+xinit;
            fy=fnc(x, params);
            if(fy==0)
                return x;
            niter++;
            if(fy*finit>0)  // go further left
                y = (y+1)/2;
        } while(fy*finit>0 && niter<60);
        if(fy*finit>0)
            return gsl_nan();
        return findRoot(fnc, params, xinit, x, reltoler);
    }
}

double findPositiveValue(function fnc, void* params, double x_0, 
                         double* out_f_0, double* out_f_p, double* out_der)
{
    double f_0 = fnc(x_0, params);
    if(out_f_0!=NULL) *out_f_0 = f_0;
    // store the initial value even if don't succeed in finding a better one
    if(out_f_p!=NULL) *out_f_p = f_0;
    if(f_0>0) {
        return x_0;
    }
    double delta = fmax(fabs(x_0) * GSL_SQRT_DBL_EPSILON, GSL_DBL_EPSILON);
    int niter=0;
    do {
        double f_minus= fnc(x_0-delta, params);
        double f_plus = fnc(x_0+delta, params);
        if(f_plus>0) {
            if(out_f_p!=NULL) *out_f_p = f_plus;
            if(out_der!=NULL) *out_der = (1.5*f_plus-2*f_0+0.5*f_minus)/delta;
            return x_0+delta;
        }
        if(f_minus>0) {
            if(out_f_p!=NULL) *out_f_p = f_minus;
            if(out_der!=NULL) *out_der = (2*f_0-1.5*f_minus-0.5*f_plus)/delta;
            return x_0-delta;
        }
        // simple recipes didn't work; estimate the first and the second derivatives..
        double der_0 = (f_plus-f_minus)/(2*delta);
        double der2_0= (f_plus+f_minus-2*f_0)/(delta*delta);
        // ..and try to guess the point of zero crossing by solving a quadratic equation
        double D=der_0*der_0-2*f_0*der2_0;  
        double x_new;
        if(D>=0) {
            D=sqrt(D);
            double x_1=x_0+(-der_0+D)/der2_0;
            double x_2=x_0+(-der_0-D)/der2_0;
            x_new = (fabs(x_1-x_0)<fabs(x_2-x_0)) ? x_1 : x_2;
        } else {  // no luck with quadratic extrapolation, try a linear one
            x_new = x_0 - (f_0-GSL_SQRT_DBL_EPSILON)/der_0;
        }
        double f_new = fnc(x_new, params);
        if(f_new>0) {
            if(out_f_p!=NULL) *out_f_p=f_new;
            // be satisfied with 1st order rule, at least the sign will be correct
            if(out_der!=NULL) *out_der = (f_new-f_0)/(x_new-x_0);
            return x_new;
        }
        if(f_new>fmax(f_minus, f_plus)) {  // move to a better location
            delta=fabs(x_0-x_new)*2;
            x_0=x_new;
            f_0=f_new;
        } else if(f_minus>f_0) { 
            x_0-=delta;
            f_0=f_minus;
        } else if(f_plus>f_0) {
            x_0+=delta;
            f_0=f_plus;
        } else
            return gsl_nan();
        niter++;
    } while(niter<3);
    return gsl_nan();
}

double integrate(function fnc, void* params, double x1, double x2, double reltoler)
{
    if(x1==x2)
        return 0;
    gsl_function F;
    F.function=fnc;
    F.params=params;
    double result, error;
    if(reltoler==0) {  // don't care about accuracy -- use the fastest integration rule
#if 1
        const int N=10;  // tables up to N=20 are hard-coded in the library, no overhead
        gsl_integration_glfixed_table* t = gsl_integration_glfixed_table_alloc(N);
        result = gsl_integration_glfixed(&F, x1, x2, t);
        gsl_integration_glfixed_table_free(t);
#else
        double dummy;  // 15-point Gauss-Kronrod
        gsl_integration_qk15(&F, x1, x2, &result, &error, &dummy, &dummy);
#endif
    } else {  // use adaptive method with limited max # of points (87)
        size_t neval;
        gsl_integration_qng(&F, x1, x2, 0, reltoler, &result, &error, &neval);
    }
    return result;
}

/** The integral \int_{x1}^{x2} f(x) dx is transformed into 
    \int_{y1}^{y2} f(x(y)) (dx/dy) dy,  where x(y) = x_low + (x_upp-x_low) y^2 (3-2y),
    and x1=x(y1), x2=x(y2).   */
struct ScaledIntParam {
    function F;
    void* param;
    double x_low, x_upp;
};

static double scaledIntegrand(double y, void* vparam) {
    ScaledIntParam* param=static_cast<ScaledIntParam*>(vparam);
    const double x = param->x_low + (param->x_upp-param->x_low) * y*y*(3-2*y);
    const double dx = (param->x_upp-param->x_low) * 6*y*(1-y);
    double val=(*(param->F))(x, param->param);
    return val*dx;
}

// invert the above relation between x and y by solving a cubic equation
static double solveForScaled_y(double x, const ScaledIntParam& param) {
    assert(x>=param.x_low && x<=param.x_upp);
    if(x==param.x_low) return 0;
    if(x==param.x_upp) return 1;
    double phi=acos(1-2*(x-param.x_low)/(param.x_upp-param.x_low))/3.0;
    return (1 - cos(phi) + M_SQRT3*sin(phi))/2.0;
}

double integrateScaled(function fnc, void* params, double x1, double x2, 
    double x_low, double x_upp, double rel_toler)
{
    if(x1==x2) return 0;
    if(x1>x2 || x1<x_low || x2>x_upp || x_low>=x_upp)
        throw std::invalid_argument("Error in integrate_scaled: arguments out of range");
    ScaledIntParam param;
    param.F=fnc;
    param.param=params;
    param.x_low=x_low;
    param.x_upp=x_upp;
    double y1=solveForScaled_y(x1, param);
    double y2=solveForScaled_y(x2, param);
    return integrate(scaledIntegrand, &param, y1, y2, rel_toler);
}

double deriv(function fnc, void* params, double x, double h, int dir)
{
    gsl_function F;
    F.function=fnc;
    F.params=params;
    double result, error;
    if(dir==0)
        gsl_deriv_central(&F, x, h, &result, &error);
    else if(dir>0)
        gsl_deriv_forward(&F, x, h, &result, &error);
    else
        gsl_deriv_backward(&F, x, h, &result, &error);
    return result;
}

double linearFitZero(unsigned int N, const double x[], const double y[])
{
    double c, cov, sumsq;
    gsl_fit_mul(x, 1, y, 1, N, &c, &cov, &sumsq);
    return c;
}

// Simple ODE integrator using Runge-Kutta Dormand-Prince 8 adaptive stepping
// dy_i/dt = f_i(t) where int (*f)(double t, const double y, double f, void *params)
struct ode_impl{
    const gsl_odeiv2_step_type * T;
    gsl_odeiv2_step * s;
    gsl_odeiv2_control * c;
    gsl_odeiv2_evolve * e;
    gsl_odeiv2_system sys;
};

OdeSolver::OdeSolver(odefunction fnc, void* params, int numvars, double abstoler, double reltoler) {
    ode_impl* data=new ode_impl;
    data->sys.function=fnc;
    data->sys.jacobian=NULL;
    data->sys.dimension=numvars;
    data->sys.params=params;
    data->s=gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk8pd, numvars);
    data->c=gsl_odeiv2_control_y_new(abstoler, reltoler);
    data->e=gsl_odeiv2_evolve_alloc(numvars);
    impl=data;
}

OdeSolver::~OdeSolver() {
    ode_impl* data=static_cast<ode_impl*>(impl);
    gsl_odeiv2_evolve_free(data->e);
    gsl_odeiv2_control_free(data->c);
    gsl_odeiv2_step_free(data->s);
    delete data;
}

int OdeSolver::advance(double tstart, double tfinish, double *y){
    ode_impl* data=static_cast<ode_impl*>(impl);
    double h=tfinish-tstart;
    double direction=(h>0?1.:-1.);
    int numstep=0;
    while ((tfinish-tstart)*direction>0 && numstep<ODE_MAX_NUM_STEP) {
        int status = gsl_odeiv2_evolve_apply (data->e, data->c, data->s, &(data->sys), &tstart, tfinish, &h, y);
        // check if computation is broken
        double test=0;
        for(unsigned int i=0; i<data->sys.dimension; i++) 
            test+=y[i];
        if (status != GSL_SUCCESS || !isFinite(test)) {
            numstep=-1;
            break;
        }
        numstep++;
    }
    if(numstep>=ODE_MAX_NUM_STEP)
        throw std::runtime_error("ODE solver: number of sub-steps exceeds maximum");
    return numstep;
}

#if 0
//=================================================================================================
// RANDOM NUMBERS //
// random number generators - rand_uniform returns random numbers uniformly distributed in the
// interval [0,1], rand_gaussian returns gaussian distributed random numbers with sd = sigma

class rand_base{
	private:
		const gsl_rng_type * TYPE;
		unsigned long int seed;
    public:
       	gsl_rng * r;
     	rand_base(unsigned long int s){
     	// construct random number generator with seed s
     		seed = s;
     		gsl_rng_env_setup();
     	    TYPE = gsl_rng_default;
       		r    = gsl_rng_alloc (TYPE);
       		gsl_rng_set(r, seed);
       		}
       	~rand_base(){gsl_rng_free (r);}
       	void reseed(unsigned long int newseed){
       	// give new seed
       		seed = newseed;
       		gsl_rng_set(r,newseed);
       		}
};

class rand_uniform:public rand_base{
	public:
		rand_uniform(unsigned long int SEED=0):rand_base(SEED){}
		~rand_uniform(){}
		// return uniformly distributed random numbers
		double nextnumber(){return gsl_rng_uniform (r);}
};

class rand_gaussian:public rand_base{
	public:
		double sigma;
		rand_gaussian(double s, unsigned long int SEED=0):rand_base(SEED){sigma = s;}
		~rand_gaussian(){}
		// return gaussian distributed random numbers
		double nextnumber(){return gsl_ran_gaussian (r,sigma);}
		void newsigma(double newsigma){sigma=newsigma;}
};

class rand_exponential:public rand_base{
    public:
        double scale;
        rand_exponential(double scale, unsigned long int SEED=0):rand_base(SEED),scale(scale){}
        ~rand_exponential(){}
        // return exponentially distributed random numbers
        double nextnumber(){return gsl_ran_exponential (r,scale);}
        void new_scale(double newscale){scale=newscale;}
};

//=================================================================================================
// ROOT FINDING	  //
// finds root by Brent's method. Constructor initialises function and tolerances, findroot finds
// root in given interval. Function must be of form double(*func)(double,void*)

class root_find{
	private:
};

//=================================================================================================
// INTEGRATION //
// Simple 1d numerical integration using adaptive Gauss-Kronrod which can deal with singularities
// constructor takes function and tolerances and integrate integrates over specified region.
// integrand function must be of the form double(*func)(double,void*)
class integrator{
	private:
		gsl_integration_workspace *w;
		void *p;
		double result,err,eps;
		gsl_function F;
		size_t neval;
	public:
		integrator(double eps): eps(eps){
			neval    = 1000;
			w        = gsl_integration_workspace_alloc (neval);
			F.params = &p;
       		}
		~integrator(){gsl_integration_workspace_free (w);}
		double integrate(double(*func)(double,void*),double xa, double xb){
		    F.function = func;
			gsl_integration_qags (&F, xa, xb, 0, eps, neval,w, &result, &err);
			//gsl_integration_qng(&F, xa, xb, 0, eps, &result, &err, &neval);
			return result;
			}
		double error(){return err;}
};

inline double integrate(double(*func)(double,void*),double xa, double xb, double eps){
	double result,err; size_t neval;void *p;gsl_function F;F.function = func;F.params = &p;
	gsl_integration_qng(&F, xa, xb, 0, eps, &result, &err, &neval);
	return result;
}

class MCintegrator{
	private:
		gsl_monte_vegas_state *s;
		const gsl_rng_type *T;
 		gsl_rng *r;
 		size_t dim;
 	public:
 		MCintegrator(size_t Dim){
 			dim=Dim;
 			gsl_rng_env_setup ();
  			T = gsl_rng_default;
  			r = gsl_rng_alloc (T);
 			s = gsl_monte_vegas_alloc(Dim);
 		}
 		~MCintegrator(){
 			gsl_monte_vegas_free(s);
 			gsl_rng_free(r);
 		}
 		double integrate(double(*func)(double*,size_t,void*),double *xlow, double *xhigh,
 		size_t calls, double *err, int burnin=10000){
 			gsl_monte_function G = { func, dim, 0 }; double res;
 			if(burnin)gsl_monte_vegas_integrate(&G,xlow,xhigh,dim,burnin,r,s,&res,err);
 			gsl_monte_vegas_integrate(&G,xlow,xhigh,dim,calls,r,s,&res,err);
 			return res;
 		}
};

//=================================================================================================
// 1D INTERPOLATION //
// Interpolation using cubic splines

class interpolator{
	private:
		gsl_interp_accel *acc;
		gsl_spline *spline;
	public:
		interpolator(double *x, double *y, int n){
			acc = gsl_interp_accel_alloc();
			spline = gsl_spline_alloc(gsl_interp_cspline,n);
			gsl_spline_init (spline, x, y, n);
		}
		~interpolator(){
			gsl_spline_free (spline);
         	gsl_interp_accel_free (acc);
        }
        double interpolate(double xi){
        	return gsl_spline_eval (spline, xi, acc);
        }
        double derivative(double xi){
        	return gsl_spline_eval_deriv(spline, xi, acc);
        }
        void new_arrays(double *x, double *y,int n){
        	spline = gsl_spline_alloc(gsl_interp_cspline,n);
        	gsl_spline_init (spline, x, y, n);
        }
};

//=================================================================================================
// SORTING //
// sorting algorithm
// sort2 sorts first argument and then applies the sorted permutation to second list
class sorter{
	private:
		const gsl_rng_type * T;
       	gsl_rng * r;
    public:
    	sorter(){
    		gsl_rng_env_setup();
            T = gsl_rng_default;
     		r = gsl_rng_alloc (T);
     		}
     	~sorter(){gsl_rng_free (r);}
     	void sort(double *data, int n){
     		gsl_sort(data,1,n);
     	}
     	void sort2(double *data, int n, double *data2){
     		size_t p[n];
     		gsl_sort_index(p,data,1,n);
     		gsl_permute(p,data2,1,n);
     	}

};

//=================================================================================================
// ODE SOLVER //


//=================================================================================================
// MINIMISER //
// finds a minimum of a function of the form double(*func)(const gsl_vector *v, void *params)
// using a downhill simplex algorithm. Setup minimiser with initial guesses and required tolerance
// with constructor and then minimise with minimise().
class minimiser{
	private:
		const gsl_multimin_fminimizer_type *T; ;
		gsl_multimin_fminimizer *s;
		gsl_vector *ss, *x;
		gsl_multimin_function minex_func;
		size_t iter; int status,N_params; double size;
		double eps;
	public:
		minimiser(double(*func)(const gsl_vector *v, void *params),double *parameters,
		 int N, double *sizes, double eps, void *params):N_params(N), eps(eps){

			T = gsl_multimin_fminimizer_nmsimplex2rand;
			ss = gsl_vector_alloc (N_params);x = gsl_vector_alloc (N_params);
			for(int i=0;i<N_params;i++){
				gsl_vector_set (x, i, parameters[i]);gsl_vector_set(ss,i,sizes[i]);}

			minex_func.n = N_params; minex_func.f = func; minex_func.params = params;
			s = gsl_multimin_fminimizer_alloc (T, N_params);
			gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
			status = 0; iter = 0;

		}

		minimiser(double(*func)(const gsl_vector *v, void *params),std::vector<double> parameters,
		 std::vector<double> sizes, double eps,void *params): eps(eps){

			N_params = parameters.size();
			T = gsl_multimin_fminimizer_nmsimplex2rand;
			ss = gsl_vector_alloc (N_params);x = gsl_vector_alloc (N_params);
			for(int i=0;i<N_params;i++){
				gsl_vector_set (x, i, parameters[i]);gsl_vector_set(ss,i,sizes[i]);}

			minex_func.n = N_params; minex_func.f = func; minex_func.params = params;
			s = gsl_multimin_fminimizer_alloc (T, N_params);
			gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
			status = 0; iter = 0;
		}

		~minimiser(){
			gsl_vector_free(x);
			gsl_vector_free(ss);
			gsl_multimin_fminimizer_free (s);
		}

		double minimise(double *results,unsigned int maxiter,bool vocal){
			do
			  {
				iter++; status = gsl_multimin_fminimizer_iterate(s);
				if(status)break;
				size = gsl_multimin_fminimizer_size (s);
				status = gsl_multimin_test_size (size, eps);
//				if(vocal){	std::cout<<iter<<" ";
//							for(int i=0; i<N_params;i++)std::cout<<gsl_vector_get(s->x,i)<<" ";
//							std::cout<<s->fval<<" "<<size<<std::endl;
//							}
			}
			while (status == GSL_CONTINUE && iter < maxiter);
			for(int i=0;i<N_params;i++){results[i] = gsl_vector_get(s->x,i);}
			return s->fval;
		}

		double minimise(std::vector<double> *results,unsigned int maxiter,bool vocal){
			do
			  {
				iter++; status = gsl_multimin_fminimizer_iterate(s);
				if(status)break;
				size = gsl_multimin_fminimizer_size (s);
				status = gsl_multimin_test_size (size, eps);
//				if(vocal){	std::cout<<iter<<" ";
//							for(int i=0; i<N_params;i++)std::cout<<gsl_vector_get(s->x,i)<<" ";
//							std::cout<<s->fval<<" "<<size<<std::endl;
//							}
			}
			while (status == GSL_CONTINUE && iter < maxiter);
			for(int i=0;i<N_params;i++) results->push_back(gsl_vector_get(s->x,i));
			return s->fval;
		}
};

class minimiser1D{
	private:
		const gsl_min_fminimizer_type *T; ;
		gsl_min_fminimizer *s;
		size_t iter; int status;
		double m, a, b, eps;
	public:
		minimiser1D(double(*func)(double, void *params), double m, double a, double b, double eps, void* params)
			:m(m), a(a), b(b), eps(eps){

			gsl_function F;F.function = func;F.params = params;
			T = gsl_min_fminimizer_brent;
			s = gsl_min_fminimizer_alloc (T);
			gsl_min_fminimizer_set (s, &F, m, a, b);
			status = 0; iter = 0;
		}
		~minimiser1D(){
			gsl_min_fminimizer_free (s);
		}
		double minimise(unsigned int maxiter){
			do
			  {
				iter++;
				status = gsl_min_fminimizer_iterate(s);
				m = gsl_min_fminimizer_x_minimum (s);
           		a = gsl_min_fminimizer_x_lower (s);
           		b = gsl_min_fminimizer_x_upper (s);
				status = gsl_min_test_interval (a, b, eps, 0.0);
			}
			while (status == GSL_CONTINUE && iter < maxiter);
			return m;
		}
};
/*
double Distance(void *xp, void *yp){
       double x = *((double *) xp);
       double y = *((double *) yp);
       return fabs(x - y);
}

void Step(const gsl_rng * r, void *xp, double step_size){
    double old_x = *((double *) xp);
    double new_x;

    double u = gsl_rng_uniform(r);
    new_x = u * 2 * step_size - step_size + old_x;

    memcpy(xp, &new_x, sizeof(new_x));
}

void Print(void *xp){
    printf ("%12g", *((double *) xp));
}

class sim_anneal{
	private:
		const gsl_rng_type * T;
    	gsl_rng * r;
    	int N_TRIES, ITER_FIXED_T;
    	double STEP_SIZE, K, T_INITIAL, MU_T, T_MIN;
    	gsl_siman_params_t params;
    public:
    	sim_anneal(int N_TRIES, int ITER_FIXED_T, double STEP_SIZE):
    	N_TRIES(N_TRIES), ITER_FIXED_T(ITER_FIXED_T),STEP_SIZE(STEP_SIZE){
    		gsl_rng_env_setup();
            T = gsl_rng_default;
     	  	r = gsl_rng_alloc(T);
     	  	//params[0]=N_TRIES;params[1]=ITER_FIZED_T;params[2]=STEP_SIZE;
     	  	//K=1.; params[3]=K; T_INITIAL=0.008; params[4]=T_INITIAL;
     	  	//MU_T=1.003; params[5]=MU_T; T_MIN=2.0e-6; params[6]=T_MIN;
     	  	gsl_siman_params_t params
       = {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
          K, T_INITIAL, MU_T, T_MIN};
    	}
    	~sim_anneal(){
    		gsl_rng_free(r);
    	}
    	double minimise(double(*func)(void *xp), double x){
    		double x_initial=x;
    		gsl_siman_solve(r, &x_initial, &func, Step, Distance, Print,
                       		NULL, NULL, NULL,
                       		sizeof(double), params);
    		return x_initial;
    	}
};*/

//=================================================================================================
// SPECIAL FUNCTIONS //
inline double erf(double x){return gsl_sf_erf (x);}
inline double erfc(double x){return gsl_sf_erfc (x);}
inline double besselI(double x, int n){return gsl_sf_bessel_In (n,x);}
inline double besselJ(double x, int n){return gsl_sf_bessel_Jn (n,x);}
inline double gamma(double x){return gsl_sf_gamma (x);}
inline double ellint_first(double phi, double k){ return gsl_sf_ellint_F(phi,k,(gsl_mode_t)1e-15);}
// F(\phi,k) = \int_0^\phi \d t \, \frac{1}{\sqrt{1-k^2\sin^2 t}}
inline double ellint_second(double phi, double k){ return gsl_sf_ellint_E(phi,k,(gsl_mode_t)1e-15);}
// E(\phi,k) = \int_0^\phi \d t \, \sqrt{1-k^2\sin^2 t}
inline double ellint_third(double phi, double k, double n){ return gsl_sf_ellint_P(phi,k,n,(gsl_mode_t)1e-15);}
// \Pi(\phi,k,n) = \int_0^\phi \d t \, \frac{1}{(1+n\sin^2 t)\sqrt{1-k^2\sin^2 t}}
#endif

}  // namespace
