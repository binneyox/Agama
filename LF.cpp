/* The next section contains classes enabling one to sample models with
 * composite DFs such that each star is assigned a magnitude as well as
 * phase-space coordinates. All classes descend from
 * DF_LFIntegrandNdim, which is slightly exteded from DFIntegrandNdim.
 * The inheriting classes are generally sipler than the classes that
 * inherit DFIntegrandNdim because they are only intended for sampling
 * rather than evaluating moments.
 */

inline double unscaleMag(const double x, const double bright, const double faint, double* jac){
	double d = (faint - bright);
	if(jac) *jac = d;
	return bright + x*d;
}

class DF_LFIntegrandNdim: public math::IFunctionNdim {
	public:
		double bright, faint;
		const obs::BaseLos* los;
		explicit DF_LFIntegrandNdim(const GalaxyModel& _model, bool separate,
					 const obs::BaseLos* _los=NULL, double _bright=NULL, double _faint=NULL)
				:
		    model(_model),
		    dflen(separate ? model.distrFunc.numValues() : 1),
		los(_los),
		bright(_bright),
		faint(_faint)
		{}

    /** compute one or more moments of distribution function. */
		virtual void eval(const double vars[], double values[]) const
		    {
	// value(s) of distribution function - more than one for a multicomponent DF treated separately
	// (array allocated on the stack, no need to delete it manually)
			    double* dfval = static_cast<double*>(alloca(dflen * sizeof(double)));
			    std::fill(dfval, dfval+dflen, 0);  // initialize with zeroes (remain so in case of errors)
	// 1. get the position/velocity components in cylindrical coordinates
			    double jac;   // jacobian of variable transformation
			    coord::PosVelCylMag posvel = unscaleVars(vars, &jac);
			    if(jac == 0) {  // we can't compute actions, but pretend that DF*jac is zero
				    outputValues(posvel, dfval, values);
				    return;
			    }

			    try{
	    // 2. determine the actions
				    actions::Actions act = model.actFinder.actions(posvel.first);

	    // 3. compute the value of distribution function times the jacobian
	    // FIXME: in some cases the Fudge action finder may fail and produce
	    // zero values of Jr,Jz instead of very large ones, which may lead to
	    // unrealistically high DF values. We therefore ignore these points
	    // entirely, but the real problem is with the action finder, not here.
				    bool valid = true;
				    if(isFinite(act.Jr + act.Jz + act.Jphi) && (act.Jr!=0 || act.Jz!=0)) {
					    if(bright){//proposed mag is apparent
						    double s = los->s(posvel.first);
						    double shift = 5*log10(100 * s);
						    shift += los->A_H(s);
						    double Mag=posvel.second - shift;
						    model.distrFunc.withLF(act, dfval, Mag);
					    } 
					    else{//proposed mag is absolute
						    model.distrFunc.withLF(act, dfval, posvel.second);  
					    }
				// multiply by jacobian, check for possibly invalid values and replace them with zeroes
					    for(unsigned int i=0; i<dflen; i++) {
						    if(!isFinite(dfval[i])) {
							    valid = false;
							    dfval[i] = 0;
						    } 
						    else
							    dfval[i] *= jac;
					    }
				    }
				    if(!valid)
					    throw std::runtime_error("DF is not finite");
			    }
			    catch(std::exception& e) {
	    // dump out the error at the highest debug logging level
				    if(utils::verbosityLevel >= utils::VL_VERBOSE) {
					    utils::msg(utils::VL_VERBOSE, "DFIntegrandNdim", std::string(e.what()) +
						    " at R="+utils::toString(posvel.first.R)  +", z="   +utils::toString(posvel.first.z)+
						    ", phi="+utils::toString(posvel.first.phi)+", vR="  +utils::toString(posvel.first.vR)+
						    ", vz=" +utils::toString(posvel.first.vz) +", vphi="+utils::toString(posvel.first.vphi));
				    }
	    // ignore the error and proceed as if the DF was zero (keep the initially assigned zero)
			    }

	// 4. output the value(s) to the integration routine
			    outputValues(posvel, dfval, values);
		    }

    /** convert from scaled variables used in the integration routine 
        to the actual position/velocity point.
        \param[in]  vars  is the array of scaled variables;
        \param[out] jac (optional)  is the jacobian of transformation, if NULL it is not computed;
        \return  the position and velocity in cylindrical coordinates.
    */
		    virtual coord::PosVelCylMag unscaleVars(const double vars[], double* jac=0) const = 0;

    /** output the value(s) computed at a given point to the integration routine.
        \param[in]  point  is the position/velocity point;
        \param[in]  dfval  is the value or array of values of distribution function at this point;
        \param[out] values is the array of one or more values that are computed
    */
		    virtual void outputValues(const coord::PosVelCylMag& point,
					      const double dfval[], double values[]) const = 0;

		    const GalaxyModel& model;  ///< reference to the galaxy model to work with
    /// number of values in the DF in the case when a composite DF has more than component, and
    /// they are requested to be considered separately, otherwise 1 (a sum of all components)
		    unsigned int dflen;
};

/* The next class is what's needed to obtain initial conditions for an
 * N-body simulation  with each star assigned an absolute magnitude.
 */
class DF_LFIntegrand6dim : public DF_LFIntegrandNdim {
	public:
		DF_LFIntegrand6dim(const GalaxyModel& _model) :
		    DF_LFIntegrandNdim(_model, false) {}

    /// input variables define 6 scaled components of position and velocity plus abs magnitude
		virtual coord::PosVelCylMag unscaleVars(const double vars[], double* jac=0) const
		{
			coord::PosVelCyl point(unscalePosVel(vars, model.potential, jac));
			const double Bright=-7, Faint=12;//brightest & faintest abs mags
			double jacMag, Mag = unscaleMag(vars[6], Bright, Faint, &jacMag);
			if(jac) (*jac) *= jacMag;
			return coord::PosVelCylMag(point, Mag);
		}

	private:
		virtual unsigned int numVars()   const { return 7; }
		virtual unsigned int numValues() const { return 1; }

    /// output contains just one value of DF (even for a multicomponent DF)
		virtual void outputValues(const coord::PosVelCylMag& ,
					  const double dfval[], double values[]) const
		{
			values[0] = dfval[0];
		}
};

/* The next class samples a LOS assigning to each star an apparent mag
 * that lies between the given limiting mags of the target survey. The
 * returned mags reflect extinction as well as distance.
 */
class DF_LFIntegrandLOS: public DF_LFIntegrandNdim {
	private:
		double vesc;                  ///< escape velocity at this position
		double zeta;                  ///< the ratio of circular to escape velocity
	public:
		DF_LFIntegrandLOS(const GalaxyModel& _model,
			       const obs::BaseLos* _los,
			       const double bright, const double faint)
				:
		    DF_LFIntegrandNdim(_model, false, _los, bright, faint)
		{
			getVesc(los->deepest(), model.potential, vesc, zeta);
		}

		    virtual coord::PosVelCylMag unscaleVars(const double vars[], double* jac=0) const
		    {
			    double jacVel, s;
			    if(los->semi_inf) s = math::unscale(math::ScalingSemiInf(), vars[0],jac);
			    else s = math::unscale(math::ScalingInf(), vars[0],jac);
			    double jacMag, mag = unscaleMag(vars[4],bright,faint,&jacMag);
			    coord::PosCyl pos = los->Rzphi(s);
			    coord::VelCyl vel = unscaleVelocity(vars+1, vesc, zeta, &jacVel);
			    if(jac) (*jac) *= jacMag * jacVel;
			    return coord::PosVelCylMag(coord::PosVelCyl(pos,vel), mag);
		    }

	// 4. output the value(s) of DF, multiplied by various combinations of velocity components:
	// {f, f*vR, f*vz, f*vphi, f*vR^2, f*vz^2, f*vphi^2, f*vR*vz, f*vR*vphi, f*vz*vphi },
	// depending on the mode of operation.
		    virtual void outputValues(const coord::PosVelCylMag& pv,
					      const double dfval[], double values[]) const
		    {
			    values[0] = dfval[0];//separate = false so only 1 value
		    }

    /// dimension of the input array (distance, 3 scaled velocity
    /// components, magnitude)
		    virtual unsigned int numVars()   const { return 5; }

    /// dimension of the output array
		    virtual unsigned int numValues() const { return 1; }
};

/* The next class enables sampling of the velocity distribution at a
 * given location and assigning to each an apparent magnitude star that
 * lies between given limiting mags (those of the target survey).
 * these magnitudes reflect extinction between the Sun and the given
 * location.
 */ 
class DF_LFIntegrandAtPoint: public DF_LFIntegrandNdim {
	private:
		const coord::PosCyl point;    ///< fixed position
		double vesc;                  ///< escape velocity at this position
		double zeta;                  ///< the ratio of circular to escape velocity
		const unsigned int numOutVal; ///< number of output values for each component of DF
	public:
		DF_LFIntegrandAtPoint(const GalaxyModel& _model,
				   const coord::PosCyl& _point,
				   obs::BaseLos* _los, double bright, double faint)
				:
		    DF_LFIntegrandNdim(_model, false, _los, bright, faint),
		    point(_point),
		    numOutVal(1)
		{
			getVesc(point, model.potential, vesc, zeta);
		}
		    
		    virtual coord::PosVelCylMag unscaleVars(const double vars[], double* jac=0) const
		    {
			    double jacMag, mag = unscaleMag(vars[3],bright,faint,&jacMag);
			    coord::PosVelCyl pv(point, unscaleVelocity(vars, vesc, zeta, jac));
			    if(jac) (*jac) *= jacMag;
			    return coord::PosVelCylMag(pv, mag);
		    }

		    virtual void outputValues(const coord::PosVelCylMag& pv,
					      const double dfval[], double values[]) const
		    {
	// 4. output the value(s) of DF, multiplied by various combinations of velocity components:
	// {f, f*vR, f*vz, f*vphi, f*vR^2, f*vz^2, f*vphi^2, f*vR*vz, f*vR*vphi, f*vz*vphi },
	// depending on the mode of operation.
			    values[0] = dfval[0];// separate = false so just one val
		    }

    /// dimension of the input array (3 scaled velocity components)
		    virtual unsigned int numVars()   const { return 4; }

    /// dimension of the output array
		    virtual unsigned int numValues() const { return 1; }
};
