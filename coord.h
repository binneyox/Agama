/** \file    coord.h 
    \brief   General-purpose coordinate types and routines
    \author  Eugene Vasiliev
    \date    2015

This module provides the general framework for working with different coordinate systems.

It is heavily templated but this shouldn't intimidate the end user, because 
the most important data structures and routines have dedicated non-templated aliases.

The fundamental data types are the following:

- coordinate systems (the simplest ones have no parameters at all);
- positions and position-velocity pairs in different coordinate systems;
- an abstract class for a scalar function defined in a particular coordinate system;
- gradients and hessians of scalar functions in different coordinate systems;
- coefficients of coordinate transformations between different systems:
  derivatives of destination coords by source coords (i.e. the jacobian matrix) and 
  second derivatives of destination coords by source coords.

The fundamental routines operating on these structures are the following:

- conversion of position and position-velocity from one coordinate system to another;
- computation of coefficients of coordinate transformation (first/second derivatives);
- transformation of gradients and hessians;
- the "all-mighty function" that uses the above primitives to perform the following task:
  suppose we have a class that computes the value, gradient and hessian of a scalar function 
  in a particular coordinate system ("evaluation CS"), and we need these quantities 
  in a different system ("output CS").
  The routine transforms the input coordinates from outputCS to evalCS, along with their 
  derivatives; computes the value, gradient and hessian in evalCS, transforms them back 
  to outputCS. A modification of this routine uses another intermediate CS for the situation 
  when a direct transformation is not implemented.
  The main application of this routine is the computation of potentials and forces 
  in different coordinate systems.
*/
#pragma once
#include "math_core.h"
#include "math_base.h"
#include "math_spline.h"

#define EXP __declspec(dllexport)

/** Classes and routines for representing position/velocity points, 
    gradients and hessians of scalar functions (e.g., gravitational potential), 
    and transformations between coordinate systems.
*/
namespace coord {

/// \name   Primitive data types: symmetry in 3d space
///@{

/** defines the symmetry properties of a function in three-dimensional space */
enum SymmetryType{ 
    ST_NONE         = 0, ///< no symmetry whatsoever
    // basic symmetries:
    ST_XREFLECTION  = 1, ///< change of sign in x (flip about yz plane)
    ST_YREFLECTION  = 2, ///< change of sign in y
    ST_ZREFLECTION  = 4, ///< change of sign in z
    ST_REFLECTION   = 8, ///< mirror reflection about origin (change of sign of all coordinates simultaneously)
    ST_ZROTATION    =16, ///< rotation about z axis
    ST_ROTATION     =32, ///< rotation about arbitrary axis
    // composite symmetries:
    /// bisymmetric: a combination of z-reflection and a mirror symmetry about origin,
    /// resulting in a mirror symmetry in xy-plane - change of sign in x and y simultaneously;
    /// this is suitable, for instance, to describe a spiral pattern containing only even-m modes
    /// (both cosine and sine terms) and symmetric w.r.t sign change in z coordinate.
    ST_BISYMMETRIC  = ST_ZREFLECTION | ST_REFLECTION,
    /// triaxial - reflection about principal planes (change of sign of any coordinate):
    /// note that while the combination of reflection symmetries about all three principal planes
    /// implies the reflection symmetry about origin (mirroring), the converse is not true, 
    /// that's why these are separate concepts; if all three plane-reflection symmetries are present,
    /// then mirror-reflection is implied, and this all is encoded in the ST_TRIAXIAL value
    ST_TRIAXIAL     = ST_XREFLECTION | ST_YREFLECTION | ST_ZREFLECTION | ST_REFLECTION, 
    ST_AXISYMMETRIC = ST_TRIAXIAL | ST_ZROTATION,    ///< axial symmetry combined with plane symmetry
    ST_SPHERICAL    = ST_AXISYMMETRIC | ST_ROTATION, ///< spherical symmetry
    ST_DEFAULT      = ST_TRIAXIAL   ///< default choice when no symmetry is specified
};

/** test for symmetry w.r.t.change of sign in x */
inline bool isXReflSymmetric(const SymmetryType sym) {
    return (sym & ST_XREFLECTION) == ST_XREFLECTION;
}

/** test for symmetry w.r.t.change of sign in y */
inline bool isYReflSymmetric(const SymmetryType sym) {
    return (sym & ST_YREFLECTION) == ST_YREFLECTION;
}

/** test for symmetry w.r.t.change of sign in z */
inline bool isZReflSymmetric(const SymmetryType sym) {
    return (sym & ST_ZREFLECTION) == ST_ZREFLECTION;
}

/** test for symmetry w.r.t.mirror reflection */
inline bool isReflSymmetric(const SymmetryType sym) {
    return (sym & ST_REFLECTION) == ST_REFLECTION;
}

/** test for rotational symmetry about z axis */
inline bool isZRotSymmetric(const SymmetryType sym) {
    return (sym & ST_ZROTATION) == ST_ZROTATION;
}

/** test for symmetry under xy-reflection */
inline bool isBisymmetric(const SymmetryType sym) {
    return (sym & ST_BISYMMETRIC) == ST_BISYMMETRIC;
}

/** test for triaxial symmetry (reflection about any of the three principal planes) */
inline bool isTriaxial(const SymmetryType sym) {
    return (sym & ST_TRIAXIAL) == ST_TRIAXIAL;
}

/** test for axisymmetry in the 'common definition'
    (i.e., invariance under rotation about z axis and under change of sign in z) */
inline bool isAxisymmetric(const SymmetryType sym) {
    return (sym & ST_AXISYMMETRIC) == ST_AXISYMMETRIC;
}

/** test for spherical symmetry */
inline bool isSpherical(const SymmetryType sym) {
    return (sym & ST_SPHERICAL) == ST_SPHERICAL;
}

///@}
/// \name   Primitive data types: coordinate systems
///@{

/// trivial coordinate systems don't have any parameters, 
/// their class names are simply used as tags in the rest of the code

/// cartesian coordinate system (galactocentric)
struct EXP Car{
    static const char* name() { return "Cartesian"; }
};

/// cylindrical coordinate system (galactocentric)
struct EXP Cyl{
    static const char* name() { return "Cylindrical"; }
};

/// spherical coordinate system (galactocentric)
struct EXP Sph{
    static const char* name() { return "Spherical"; }
};

/// spherical coordinate system with modified polar angle variable
struct EXP SphMod{
    static const char* name() { return "Modified spherical"; }
};

//  less trivial:
/** prolate spheroidal coordinate system, defined by a single parameter 
    Delta>0 (focal distance).
    The traditionally used two parameters alpha and gamma (e.g., de Zeeuw 1985) 
    are not independent, so we define  Delta^2 = gamma - alpha.
*/
struct EXP ProlSph{
	const double Delta2;     ///< Delta^2 = gamma - alpha > 0
	ProlSph(double Delta);   ///< Delta is the focal distance
	static const char* name() { return "Prolate spheroidal"; }
};

struct EXP ProlMod{
    const double D;
    ProlMod(double _D) : D(_D) {};
    static const char* name() { return "Modified prolate spheroidal"; }
};

struct EXP UVSph{
	double Delta, Delta2;     ///< Delta^2 = gamma - alpha > 0
	UVSph(double D){
		if(Delta<=0)
			throw std::invalid_argument("Invalid parameters for UV Spheroidal coordinate system");
		Delta = D; Delta2 = D*D;
	}
	UVSph(){}
	void set(double D){
		Delta=D; Delta2=D*D;
	}
	static const char* name() { return "UV Prolate spheroidal"; }
};


///@}
/// \name   Primitive data types: position in different coordinate systems
///@{

/// position in arbitrary coordinates:
/// the data types are defined as templates with the template parameter
/// being any of the coordinate system names defined above
template<typename CoordT> struct EXP PosT;

/// position in cartesian coordinates
template<> struct EXP PosT<Car>{
    double x, y, z;   ///< three cartesian coordinates
    PosT<Car>() {};
    PosT<Car>(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {};
};
/// an alias to templated type specialization of position in cartesian coordinates
typedef struct EXP PosT<Car> PosCar;


/// position in cylindrical coordinates
template<> struct EXP PosT<Cyl>{
    double R;   ///< cylindrical radius = sqrt(x^2+y^2)
    double z;   ///< z coordinate
    double phi; ///< azimuthal angle in x-y plane [0:2pi)
    PosT<Cyl>() {};
    PosT<Cyl>(double _R, double _z, double _phi) : R(_R), z(_z), phi(_phi) {};
};
typedef struct EXP PosT<Cyl> PosCyl;

/// position in spherical coordinates
template<> struct EXP PosT<Sph>{
    double r;     ///< spherical radius
    double theta; ///< polar angle [0:pi] - 0 means along z axis in positive direction,
                  ///< pi is along z in negative direction, pi/2 is in x-y plane
    double phi;   ///< azimuthal angle in x-y plane [0:2pi)
    PosT<Sph>() {};
    PosT<Sph>(double _r, double _theta, double _phi) : r(_r), theta(_theta), phi(_phi) {};
};
typedef struct EXP PosT<Sph> PosSph;

/// position in modified spherical coordinates
template<> struct PosT<SphMod>{
    double r;   ///< spherical radius
    double tau; ///< replacement for polar angle theta: tau = cos(theta) / (1 + sin(theta)); -1<=tau<=1
    double phi; ///< azimuthal angle in x-y plane [0:2pi)
    PosT<SphMod>() {};
    PosT<SphMod>(double _r, double _tau, double _phi) : r(_r), tau(_tau), phi(_phi) {};
};
typedef struct EXP PosT<SphMod> PosSphMod;

/** position in prolate spheroidal coordinates.
    We use a somewhat different definition from de Zeeuw 1985, namely: 
    the value of `nu` keeps track of the sign of z, so that the conversion between cylindrical
    and prolate spheroidal coordinates is invertible. */
template<> struct EXP PosT<ProlSph>{
	double lambda;  ///< lies in the range [delta:infinity)
	double nu;      ///< lies in the range [-delta:delta]; negative for z<0
	double phi;     ///< usual azimuthal angle
	const ProlSph& coordsys;  ///< a point means nothing without specifying its coordinate system
	PosT<ProlSph>(double _lambda, double _nu, double _phi, const ProlSph& _coordsys):
	    lambda(_lambda), nu(_nu), phi(_phi), coordsys(_coordsys) {};
};
typedef struct EXP PosT<ProlSph> PosProlSph;

/** position in prolate spheroidal coordinates (alternative version).
    The original formulation in terms of u and v variables, such that the cylindrical coordinates
    are given by  R = D sinh(u) sin(v), z = D cosh(u) cos(v),
    is replaced by  rho = D sinh(u), tau = cos(v) / (1 + sin(v)).
    Here D is the focal distance of the coordinate system, which may even be zero.
    An auxiliary variable chi is set equal to D cosh(u) = sqrt(D^2 + rho^2);
    as all member variables are declared const, this assignment may not be accidentally changed.
*/
template<> struct EXP PosT<ProlMod>{
    const double rho;  ///< lies in the range [0:infinity), equal to cylindrical radius when z=0
    const double tau;  ///< lies in the range [-1:1], analog of tau in SphMod coords
    const double phi;  ///< usual azimuthal angle
    const double chi;  ///< equal to sqrt(rho^2 + D^2), where D is the focal distance of coord.sys.
    PosT<ProlMod>(double _rho, double _tau, double _phi) : 
        rho(_rho), tau(_tau), phi(_phi), chi(_rho) {};
    PosT<ProlMod>(double _rho, double _tau, double _phi, double _chi) : 
        rho(_rho), tau(_tau), phi(_phi), chi(_chi) {};
    PosT<ProlMod>(double _rho, double _tau, double _phi, const ProlMod& coordsys);
};
typedef struct EXP PosT<ProlMod> PosProlMod;

template<> struct EXP PosT<UVSph>{
	double u;  ///< lies in the range [0:infinity)
	double v;      ///< lies in the range [0,pi]; z=0 <-> v=pi/2
	double phi;     ///< usual azimuthal angle
	const UVSph& coordsys;  ///< a point means nothing without specifying its coordinate system
	PosT<UVSph>(double _u, double _v, double _phi, const UVSph& _coordsys):
	    u(_u), v(_v), phi(_phi), coordsys(_coordsys) {};
};
typedef struct EXP PosT<UVSph> PosUVSph;
///@}
/// \name   Primitive data types: velocity in different coordinate systems
///@{

/// velocity in arbitrary coordinates
template<typename CoordT> struct EXP VelT;

/// velocity in cartesian coordinates
template<> struct EXP VelT<Car> {
	double vx, vy, vz;   ///< components of velocity along three cartesian axes
	VelT<Car>() {};
	VelT<Car>(double _vx, double _vy, double _vz) : vx(_vx), vy(_vy), vz(_vz) {};
};
/// an alias to templated type specialization of velocity for cartesian coordinates
typedef struct EXP VelT<Car> VelCar;

/// velocity in cylindrical coordinates 
/// (this is not the same as time derivative of position in these coordinates!)
template<> struct EXP VelT<Cyl> {
	double vR, vz, vphi;
	VelT<Cyl>() {};
	VelT<Cyl>(double _vR, double _vz, double _vphi) : vR(_vR), vz(_vz), vphi(_vphi) {};
};
typedef struct EXP VelT<Cyl> VelCyl;

/// velocity in spherical coordinates
/// (this is not the same as time derivative of position in these coordinates!)
template<> struct EXP VelT<Sph> {
	double vr, vtheta, vphi;
	VelT<Sph>() {};
	VelT<Sph>(double _vr, double _vtheta, double _vphi) : vr(_vr), vtheta(_vtheta), vphi(_vphi) {};
};
typedef struct EXP VelT<Sph> VelSph;

/// momentum in arbitrary coordinates
template<typename CoordT> struct EXP MomT;

/// moentum in cartesian coordinates
template<> struct EXP MomT<Car> {
	double px, py, pz;   ///< components of momentum along three cartesian axes
	MomT<Car>() {};
	MomT<Car>(double _px, double _py, double _pz) : px(_px), py(_py), pz(_pz) {};
};
/// an alias to templated type specialization of momentum for cartesian coordinates
typedef struct EXP MomT<Car> MomCar;

/// momentum in cylindrical coordinates 
template<> struct EXP MomT<Cyl> {
	double pR, pz, pphi;
	MomT<Cyl>() {};
	MomT<Cyl>(double _pR, double _pz, double _pphi) : pR(_pR), pz(_pz), pphi(_pphi) {};
};
typedef struct EXP MomT<Cyl> MomCyl;

/// momentum in spherical coordinates
template<> struct EXP MomT<Sph> {
	double pr, ptheta, pphi;
	MomT<Sph>() {};
	MomT<Sph>(double _pr, double _ptheta, double _pphi) : pr(_pr), ptheta(_ptheta), pphi(_pphi) {};
};
typedef struct EXP MomT<Sph> MomSph;

/// momentum in prolate spheroidal coordinates, canonically conjugate to the position
template<> struct EXP MomT<ProlMod> {
	double prho, ptau, pphi;
	MomT<ProlMod>() {};
	MomT<ProlMod>(double _prho, double _ptau, double _pphi) : prho(_prho), ptau(_ptau), pphi(_pphi) {};
};
typedef struct EXP MomT<ProlMod> MomProlMod;

/// momentum in UV prolate spheroidal coordinates, canonically conjugate to the position
template<> struct EXP MomT<UVSph> {
	double pu, pv, pphi;
	MomT<UVSph>() {};
	MomT<UVSph>(double _pu, double _pv, double _pphi) : pu(_pu), pv(_pv), pphi(_pphi) {};
};
typedef struct EXP MomT<UVSph> MomUVSph;

///@}
/// \name   Primitive data types: second moments of velocity in different coordinate systems
///@{

/// second moment of velocity in arbitrary coordinates
template<typename CoordT> struct EXP Vel2T;

/// velocity in cartesian coordinates
template<> struct EXP Vel2T<Car> {
    double vx2, vy2, vz2, vxvy, vxvz, vyvz;
};
typedef struct EXP Vel2T<Car> Vel2Car;

/// second moment of velocity in cylindrical coordinates 
template<> struct EXP Vel2T<Cyl> {
    double vR2, vz2, vphi2, vRvz, vRvphi, vzvphi;
};
typedef struct EXP Vel2T<Cyl> Vel2Cyl;

/// second moment of velocity in spherical coordinates
template<> struct EXP Vel2T<Sph> {
    double vr2, vtheta2, vphi2, vrvtheta, vrvphi, vthetavphi;
};
typedef struct EXP Vel2T<Sph> Vel2Sph;

///@}
/// \name   Primitive data types: position-velocity pairs in different coordinate systems
///@{

/// combined position and velocity in arbitrary coordinates
template<typename CoordT> struct EXP PosVelT;

/// combined position and velocity in cartesian coordinates
template<> struct EXP PosVelT<Car>: public PosCar, public VelCar {
	PosVelT<Car>() {};
    /// initialize from position and velocity
	PosVelT<Car>(const PosCar& pos, const VelCar& vel) : PosCar(pos), VelCar(vel) {}
    /// initialize from explicitly given numbers
	PosVelT<Car>(double _x, double _y, double _z, double _vx, double _vy, double _vz) :
	    PosCar(_x, _y, _z), VelCar(_vx, _vy, _vz) {}
    /// initialize from an array of 6 floats (i.e., from a serialized array)
	PosVelT<Car>(const double p[]) :
	    PosCar(p[0], p[1], p[2]), VelCar(p[3], p[4], p[5]) {}
    /// serialize into an array of 6 floating-point numbers
	void unpack_to(double *out) const {
		out[0]=x; out[1]=y; out[2]=z; out[3]=vx; out[4]=vy; out[5]=vz; }
};
/// an alias to templated type specialization of position and velocity for cartesian coordinates
typedef struct EXP PosVelT<Car> PosVelCar;

/// combined position and velocity in cylindrical coordinates
template<> struct EXP PosVelT<Cyl>: public PosCyl, public VelCyl {
	PosVelT<Cyl>() {};
    /// initialize from position and velocity
	PosVelT<Cyl>(const PosCyl& pos, const VelCyl& vel) : PosCyl(pos), VelCyl(vel) {}
    /// initialize from explicitly given numbers
	PosVelT<Cyl>(double _R, double _z, double _phi, double _vR, double _vz, double _vphi) :
	    PosCyl(_R, _z, _phi), VelCyl(_vR, _vz, _vphi) {};
    /// initialize from an array of 6 floats (i.e., from a serialized array)
	PosVelT<Cyl>(const double p[]) :
	    PosCyl(p[0], p[1], p[2]), VelCyl(p[3], p[4], p[5]) {};
    /// serialize into an array of 6 floating-point numbers
	void unpack_to(double *out) const {
		out[0]=R; out[1]=z; out[2]=phi; out[3]=vR; out[4]=vz; out[5]=vphi; }
};
typedef struct EXP PosVelT<Cyl> PosVelCyl;
typedef std::pair<VelCyl, double> VelCylMag;
typedef std::pair<PosVelCyl, double> PosVelCylMag;

/// combined position and velocity in spherical coordinates
template<> struct EXP PosVelT<Sph>: public PosSph, public VelSph {
	PosVelT<Sph>() {};
    /// initialize from position and velocity
	PosVelT<Sph>(const PosSph& pos, const VelSph& vel) : PosSph(pos), VelSph(vel) {}
    /// initialize from explicitly given numbers
	PosVelT<Sph>(double _r, double _theta, double _phi, double _vr, double _vtheta, double _vphi) :
	    PosSph(_r, _theta, _phi), VelSph(_vr, _vtheta, _vphi) {};
    /// initialize from an array of 6 floats (i.e., from a serialized array)
	PosVelT<Sph>(const double p[]) :
	    PosSph(p[0], p[1], p[2]), VelSph(p[3], p[4], p[5]) {};
    /// serialize into an array of 6 floating-point numbers
	void unpack_to(double *out) const {
		out[0]=r; out[1]=theta; out[2]=phi; out[3]=vr; out[4]=vtheta; out[5]=vphi; }
};
typedef struct EXP PosVelT<Sph> PosVelSph;

/// position and velocity in prolate spheroidal coordinates
template<> struct EXP PosVelT<ProlSph>: public PosProlSph{
	double lambdadot, nudot, phidot;  ///< time derivatives of position variables
	PosVelT<ProlSph>(const PosProlSph& pos, double _lambdadot, double _nudot, double _phidot):
	    PosProlSph(pos), lambdadot(_lambdadot), nudot(_nudot), phidot(_phidot) {};
	void unpack_to(double *out) const {
		out[0]=lambda; out[1]=nu; out[2]=phi; out[3]=lambdadot; out[4]=nudot; out[5]=phidot; }
};
typedef struct EXP PosVelT<ProlSph> PosVelProlSph;

/// combined position and momentum in arbitrary coordinates
template<typename CoordT> struct EXP PosMomT;

/// canonically conjugate coordinate and momenta in modified prolate spherical coordinates
template<> struct EXP PosMomT<ProlMod>: public PosProlMod, public MomProlMod {
    PosMomT<ProlMod>(const PosProlMod& pos, const MomProlMod& mom) : PosProlMod(pos), MomProlMod(mom) {}
};
typedef struct EXP PosMomT<ProlMod> PosMomProlMod;

/// combined position and momentum in cartesian coordinates
template<> struct EXP PosMomT<Car>: public PosCar, public MomCar {
	PosMomT<Car>() {};
    /// initialize from position and velocity
	PosMomT<Car>(const PosCar& pos, const MomCar& mom) : PosCar(pos), MomCar(mom) {}
    /// initialize from explicitly given numbers
	PosMomT<Car>(double _x, double _y, double _z, double _px, double _py, double _pz) :
	    PosCar(_x, _y, _z), MomCar(_px, _py, _pz) {}
    /// initialize from an array of 6 floats (i.e., from a serialized array)
	PosMomT<Car>(const double p[]) :
	    PosCar(p[0], p[1], p[2]), MomCar(p[3], p[4], p[5]) {}
    /// serialize into an array of 6 floating-point numbers
	void unpack_to(double *out) const {
		out[0]=x; out[1]=y; out[2]=z; out[3]=px; out[4]=py; out[5]=pz; }
};
/// an alias to templated type specialization of position and velocity for cartesian coordinates
typedef struct EXP PosMomT<Car> PosMomCar;

/// combined position and momentum in cylindrical coordinates
template<> struct EXP PosMomT<Cyl>: public PosCyl, public MomCyl {
	PosMomT<Cyl>() {};
    /// initialize from position and velocity
	PosMomT<Cyl>(const PosCyl& pos, const MomCyl& vel) : PosCyl(pos), MomCyl(vel) {}
    /// initialize from explicitly given numbers
	PosMomT<Cyl>(double _R, double _z, double _phi, double _pR, double _pz, double _pphi) :
	    PosCyl(_R, _z, _phi), MomCyl(_pR, _pz, _pphi) {};
    /// initialize from an array of 6 floats (i.e., from a serialized array)
	PosMomT<Cyl>(const double p[]) :
	    PosCyl(p[0], p[1], p[2]), MomCyl(p[3], p[4], p[5]) {};
    /// serialize into an array of 6 floating-point numbers
	void unpack_to(double *out) const {
		out[0]=R; out[1]=z; out[2]=phi; out[3]=pR; out[4]=pz; out[5]=pphi; }
};
typedef struct EXP PosMomT<Cyl> PosMomCyl;

/// combined position and momentum in spherical coordinates
template<> struct EXP PosMomT<Sph>: public PosSph, public MomSph {
	PosMomT<Sph>() {};
    /// initialize from position and velocity
	PosMomT<Sph>(const PosSph& pos, const MomSph& mom) : PosSph(pos), MomSph(mom) {}
    /// initialize from explicitly given numbers
	PosMomT<Sph>(double _r, double _theta, double _phi, double _vr, double _vtheta, double _vphi) :
	    PosSph(_r, _theta, _phi), MomSph(_vr, _vtheta, _vphi) {};
    /// initialize from an array of 6 floats (i.e., from a serialized array)
	PosMomT<Sph>(const double p[]) :
	    PosSph(p[0], p[1], p[2]), MomSph(p[3], p[4], p[5]) {};
    /// serialize into an array of 6 floating-point numbers
	void unpack_to(double *out) const {
		out[0]=r; out[1]=theta; out[2]=phi; out[3]=pr; out[4]=ptheta; out[5]=pphi; }
};
typedef struct EXP PosMomT<Sph> PosMomSph;

/// canonically conjugate coordinate and momenta in modified spherical coordinates
template<> struct EXP PosVelT<SphMod>: public PosSphMod {
	double pr;   ///< p_r   = v_r = dr/dt
	double ptau; ///< p_tau = -2 * r * v_theta / (1+tau^2)
	double pphi; ///< p_phi = R * v_phi
	PosVelT<SphMod>() {};
	PosVelT<SphMod>(double _r, double _tau, double _phi, double _pr, double _ptau, double _pphi) :
	    PosSphMod(_r, _tau, _phi), pr(_pr), ptau(_ptau), pphi(_pphi) {};
};
typedef struct EXP PosVelT<SphMod> PosVelSphMod;

/// position and velocity in UV prolate spheroidal coordinates
template<> struct EXP PosVelT<UVSph>: public PosUVSph{
	double udot, vdot, phidot;  ///< time derivatives of position variables
	PosVelT<UVSph>(const PosUVSph& pos, double _udot, double _vdot, double _phidot):
	    PosUVSph(pos), udot(_udot), vdot(_vdot), phidot(_phidot) {};
	void unpack_to(double *out) const {
		out[0]=u; out[1]=v; out[2]=phi; out[3]=udot; out[4]=vdot; out[5]=phidot; }
};
typedef struct EXP PosVelT<UVSph> PosVelUVSph;

/// position and velocity in UV prolate spheroidal coordinates
template<> struct EXP PosMomT<UVSph>: public PosUVSph, public MomUVSph {
	double udot, vdot, phidot;  ///< time derivatives of position variables
	PosMomT<UVSph>(const PosUVSph& pos, const MomUVSph& mom):
	    PosUVSph(pos), MomUVSph(mom) {};
	void unpack_to(double *out) const {
		out[0]=u; out[1]=v; out[2]=phi; out[3]=pu; out[4]=pv; out[5]=pphi; }
};
typedef struct EXP PosMomT<UVSph> PosMomUVSph;

///@}
/// \name   Primitive data types: gradient of a scalar function in different coordinate systems
///@{

/// components of a gradient in a given coordinate system
template<typename CoordT> struct EXP GradT;

/// gradient of scalar function in cartesian coordinates
template<> struct EXP GradT<Car>{
    double dx, dy, dz;
};
/// an alias to templated type specialization of gradient for cartesian coordinates
typedef struct EXP GradT<Car> GradCar;

/// gradient of scalar function in cylindrical coordinates
template<> struct EXP GradT<Cyl>{
    double dR, dz, dphi;
};
typedef struct EXP GradT<Cyl> GradCyl;

/// gradient of scalar function in spherical coordinates
template<> struct EXP GradT<Sph>{
    double dr, dtheta, dphi;
};
typedef struct EXP GradT<Sph> GradSph;

/// gradient of scalar function in prolate spheroidal coordinates
template<> struct EXP GradT<ProlSph>{
	double dlambda, dnu, dphi;
};
typedef struct EXP GradT<ProlSph> GradProlSph;

template<> struct EXP GradT<ProlMod>{
    double drho, dtau, dphi;
};
typedef struct EXP GradT<ProlMod> GradProlMod;

template<> struct EXP GradT<UVSph>{
	double du, dv, dphi;
};
typedef struct EXP GradT<UVSph> GradUVSph;

///@}
/// \name   Primitive data types: hessian of a scalar function in different coordinate systems
///@{

/// components of a hessian of a scalar function (matrix of its second derivatives)
template<typename CoordT> struct EXP HessT;

/// Hessian of scalar function F in cartesian coordinates: d2F/dx^2, d2F/dxdy, etc
template<> struct EXP HessT<Car>{
    double dx2, dy2, dz2, dxdy, dydz, dxdz;
};
typedef struct EXP HessT<Car> HessCar;

/// Hessian of scalar function in cylindrical coordinates
template<> struct EXP HessT<Cyl>{
    double dR2, dz2, dphi2, dRdz, dzdphi, dRdphi;
};
typedef struct EXP HessT<Cyl> HessCyl;

/// Hessian of scalar function in spherical coordinates
template<> struct EXP HessT<Sph>{
    double dr2, dtheta2, dphi2, drdtheta, dthetadphi, drdphi;
};
typedef struct EXP HessT<Sph> HessSph;

/// Hessian of scalar function in prolate spheroidal coordinates
template<> struct EXP HessT<ProlSph>{
	double dlambda2, dnu2, dlambdadnu;  ///< note: derivatives by phi are assumed to be zero
};
typedef struct EXP HessT<ProlSph> HessProlSph;

/// Hessian of scalar function in modified prolate spheroidal coordinates
template<> struct EXP HessT<ProlMod>{
    double drho2, dtau2, drhodtau;  ///< note: derivatives by phi are assumed to be zero
};
typedef struct EXP HessT<ProlMod> HessProlMod;

template<> struct EXP HessT<UVSph>{
	double du2, dv2, dudv;  ///< note: derivatives by phi are assumed to be zero
};
typedef struct EXP HessT<UVSph> HessUVSph;


///@}
/// \name   Abstract interface classes for scalar functions
///@{

/** Prototype of a scalar function which is computed in a particular coordinate system */
template<typename CoordT>
class EXP IScalarFunction {
public:
    IScalarFunction() {};
    virtual ~IScalarFunction() {};
    /** Evaluate any combination of value, gradient and hessian of the function at a given point.
        Each of these quantities is computed and stored in the output pointer if it was not NULL. */
    virtual void evalScalar(const PosT<CoordT>& x,
        double* value=NULL,
        GradT<CoordT>* deriv=NULL,
        HessT<CoordT>* deriv2=NULL) const=0;
};

///@}
/// \name   Data types containing conversion coefficients between different coordinate systems
///@{

/** derivatives of coordinate transformation from source to destination 
    coordinate systems (srcCS=>destCS): derivatives of destination variables 
    w.r.t.source variables, aka Jacobian */
template<typename srcCS, typename destCS> struct EXP PosDerivT;

/** instantiations of the general template for derivatives of coordinate transformations
    are separate structures for each pair of coordinate systems */
template<> struct EXP PosDerivT<Car, Cyl> {
    double dRdx, dRdy, dphidx, dphidy;
};
template<> struct EXP PosDerivT<Car, Sph> {
    double drdx, drdy, drdz, dthetadx, dthetady, dthetadz, dphidx, dphidy;
};
template<> struct EXP PosDerivT<Cyl, Car> {
    double dxdR, dxdphi, dydR, dydphi;
};
template<> struct EXP PosDerivT<Cyl, Sph> {
    double drdR, drdz, dthetadR, dthetadz;
};
template<> struct EXP PosDerivT<Sph, Car> {
    double dxdr, dxdtheta, dxdphi, dydr, dydtheta, dydphi, dzdr, dzdtheta;
};
template<> struct EXP PosDerivT<Sph, Cyl> {
    double dRdr, dRdtheta, dzdr, dzdtheta;
};
template<> struct EXP PosDerivT<Cyl, ProlSph> {
	double dlambdadR, dlambdadz, dnudR, dnudz;
};
template<> struct EXP PosDerivT<ProlSph, Cyl> {
	double dRdlambda, dRdnu, dzdlambda, dzdnu;
};
template<> struct EXP PosDerivT<Cyl, ProlMod> {
    double drhodR, drhodz, dtaudR, dtaudz;
};
template<> struct EXP PosDerivT<ProlMod, Cyl> {
    double dRdrho, dRdtau, dzdrho, dzdtau;
};
template<> struct EXP PosDerivT<Cyl, UVSph> {
	double dudR, dudz, dvdR, dvdz;
};
template<> struct EXP PosDerivT<UVSph, Cyl> {
	double dRdu, dRdv, dzdu, dzdv;
};


/** second derivatives of coordinate transformation from source to destination 
    coordinate systems (srcCS=>destCS): d^2(dest_coord)/d(source_coord1)d(source_coord2) */
template<typename srcCS, typename destCS> struct EXP PosDeriv2T;

/** instantiations of the general template for second derivatives of coordinate transformations */
template<> struct EXP PosDeriv2T<Cyl, Car> {
    double d2xdRdphi, d2xdphi2, d2ydRdphi, d2ydphi2;
};
template<> struct EXP PosDeriv2T<Sph, Car> {
    double d2xdrdtheta, d2xdrdphi, d2xdtheta2, d2xdthetadphi, d2xdphi2,
           d2ydrdtheta, d2ydrdphi, d2ydtheta2, d2ydthetadphi, d2ydphi2,
           d2zdrdtheta, d2zdtheta2;
};
template<> struct EXP PosDeriv2T<Car, Cyl> {
    double d2Rdx2, d2Rdxdy, d2Rdy2, d2phidx2, d2phidxdy, d2phidy2;
};
template<> struct EXP PosDeriv2T<Sph, Cyl> {
    double d2Rdrdtheta, d2Rdtheta2, d2zdrdtheta, d2zdtheta2;
};
template<> struct EXP PosDeriv2T<Car, Sph> {
    double d2rdx2, d2rdxdy, d2rdxdz, d2rdy2, d2rdydz, d2rdz2,
        d2thetadx2, d2thetadxdy, d2thetadxdz, d2thetady2, d2thetadydz, d2thetadz2,
        d2phidx2, d2phidxdy, d2phidy2;
};
template<> struct EXP PosDeriv2T<Cyl, Sph> {
    double d2rdR2, d2rdRdz, d2rdz2, d2thetadR2, d2thetadRdz, d2thetadz2;
};
template<> struct EXP PosDeriv2T<Cyl, ProlSph> {
	double d2lambdadR2, d2lambdadRdz, d2lambdadz2, d2nudR2, d2nudRdz, d2nudz2;
};
template<> struct EXP PosDeriv2T<ProlSph, Cyl> {
	double d2Rdlambda2, d2Rdlambdadnu, d2Rdnu2, d2zdlambda2, d2zdlambdadnu, d2zdnu2;
};
template<> struct EXP PosDeriv2T<ProlMod, Cyl>{
    double d2Rdrho2, d2Rdrhodtau, d2Rdtau2, d2zdrho2, d2zdrhodtau, d2zdtau2;
};
template<> struct EXP PosDeriv2T<Cyl, ProlMod>{
    double d2rhodR2, d2rhodRdz, d2rhodz2, d2taudR2, d2taudRdz, d2taudz2;
};
template<> struct EXP PosDeriv2T<Cyl, UVSph> {
	double d2udR2, d2udRdz, d2udz2, d2vdR2, d2vdRdz, d2vdz2;
};
template<> struct EXP PosDeriv2T<UVSph, Cyl> {
	double d2Rdu2, d2Rdudv, d2Rdv2, d2zdu2, d2zdudv, d2zdv2;
};

///@}
/// \name   Routines for conversion between position/velocity in different coordinate systems
///@{

/** universal templated conversion function for positions:
    template parameters srcCS and destCS may be any of the coordinate system names.
    This template function shouldn't be used directly, because the return type depends 
    on the template and hence cannot be automatically inferred by the compiler. 
    Instead, named functions for each target coordinate system are defined below. */
template<typename srcCS, typename destCS>
EXP PosT<destCS> toPos(const PosT<srcCS>& from);

/** templated conversion taking the parameters of coordinate system into account */
template<typename srcCS, typename destCS>
PosT<destCS> toPos(const PosT<srcCS>& from, const destCS& coordsys);

/** templated conversion functions for positions 
    with names reflecting the target coordinate system. */
template<typename srcCS>
inline PosCar toPosCar(const PosT<srcCS>& from) { return toPos<srcCS, Car>(from); }
template<typename srcCS>
inline PosCyl toPosCyl(const PosT<srcCS>& from) { return toPos<srcCS, Cyl>(from); }
template<typename srcCS>
inline PosSph toPosSph(const PosT<srcCS>& from) { return toPos<srcCS, Sph>(from); }
template<typename srcCS>//dont we need  here?
inline PosUVSph toPosUVSph(const PosT<srcCS>& from) { return toPos<srcCS, UVSph>(from); }



/** universal templated conversion function for coordinates and velocities:
    template parameters srcCS and destCS may be any of the coordinate system names */
template<typename srcCS, typename destCS>
EXP PosVelT<destCS> toPosVel(const PosVelT<srcCS>& from);

/** templated conversion functions for coordinates and velocities
    with names reflecting the target coordinate system. */
template<typename srcCS>
inline PosVelCar toPosVelCar(const PosVelT<srcCS>& from) { return toPosVel<srcCS, Car>(from); }
template<typename srcCS>
inline PosVelCyl toPosVelCyl(const PosVelT<srcCS>& from) { return toPosVel<srcCS, Cyl>(from); }
template<typename srcCS>
inline PosVelSph toPosVelSph(const PosVelT<srcCS>& from) { return toPosVel<srcCS, Sph>(from); }

/** templated conversion taking the parameters of coordinate system into account */
template<typename srcCS, typename destCS>
EXP PosVelT<destCS> toPosVel(const PosVelT<srcCS>& from, const destCS& coordsys);

/** trivial conversions */
template<> inline PosCar toPos<Car,Car>(const PosCar& p) { return p;}
template<> inline PosCyl toPos<Cyl,Cyl>(const PosCyl& p) { return p;}
template<> inline PosSph toPos<Sph,Sph>(const PosSph& p) { return p;}
template<> inline PosVelCar toPosVel<Car,Car>(const PosVelCar& p) { return p;}
template<> inline PosVelCyl toPosVel<Cyl,Cyl>(const PosVelCyl& p) { return p;}
template<> inline PosVelSph toPosVel<Sph,Sph>(const PosVelSph& p) { return p;}

template<typename srcCS, typename destCS>
EXP PosMomT<destCS> toPosMom(const PosMomT<srcCS>& from);

/** templated conversion functions for coordinates and velocities
    with names reflecting the target coordinate system. */
template<typename srcCS>
inline PosMomCar toPosMomCar(const PosMomT<srcCS>& from) { return toPosMom<srcCS, Car>(from); }
template<typename srcCS>
inline PosMomCyl toPosMomCyl(const PosMomT<srcCS>& from) { return toPosMom<srcCS, Cyl>(from); }
template<typename srcCS>
inline PosMomSph toPosMomSph(const PosMomT<srcCS>& from) { return toPosMom<srcCS, Sph>(from); }

/** templated conversion taking the parameters of coordinate system into account */
template<typename srcCS, typename destCS>
EXP PosMomT<destCS> toPosMom(const PosMomT<srcCS>& from, const destCS& coordsys);

/** trivial conversions */
template<> inline PosMomCar toPosMom<Car,Car>(const PosMomCar& p) { return p;}
template<> inline PosMomCyl toPosMom<Cyl,Cyl>(const PosMomCyl& p) { return p;}
template<> inline PosMomSph toPosMom<Sph,Sph>(const PosMomSph& p) { return p;}

/** more templated conversion taking the parameters of coordinate
 ** system into account and involving both Vel and Mom*/
template<typename srcCS, typename destCS>
EXP PosMomT<destCS> toPosMom(const PosVelT<srcCS>& from, const destCS& coordsys);
template<typename srcCS, typename destCS>
EXP PosVelT<destCS> toPosVel(const PosMomT<srcCS>& from);

EXP void UVmomenta(const PosUVSph&,const VelCyl&,double&,double&);
EXP void Rzmomenta(const PosUVSph&,const double,const double,VelCyl&);


///@}
/// \name   Routines for conversion between position in different coordinate systems with derivatives
///@{

/** universal templated function for coordinate conversion that provides derivatives of transformation.
    Template parameters srcCS and destCS may be any of the coordinate system names;
    \param[in]  from specifies the point in srcCS coordinate system;
    \param[out] deriv will contain derivatives of the transformation 
                (destination coords over source coords);
    \param[out] deriv2 if not NULL, will contain second derivatives of the coordinate transformation;
    \return     point in destCS coordinate system. */
template<typename srcCS, typename destCS>
EXP PosT<destCS> toPosDeriv(const PosT<srcCS>& from, 
    PosDerivT<srcCS, destCS>* deriv, PosDeriv2T<srcCS, destCS>* deriv2=NULL);

/** templated conversion with derivatives, taking the parameters of coordinate system into account */
template<typename srcCS, typename destCS>
EXP PosT<destCS> toPosDeriv(const PosT<srcCS>& from, const destCS& coordsys,
    PosDerivT<srcCS, destCS>* deriv, PosDeriv2T<srcCS, destCS>* deriv2=NULL);

///@}
/// \name   Routines for conversion of gradients and hessians between coordinate systems
///@{

/** templated function for transforming a gradient to a different coordinate system */
template<typename srcCS, typename destCS>
EXP GradT<destCS> toGrad(const GradT<srcCS>& src, const PosDerivT<destCS, srcCS>& deriv);

/** templated function for transforming a hessian to a different coordinate system */
template<typename srcCS, typename destCS>
EXP HessT<destCS> toHess(const GradT<srcCS>& srcGrad, const HessT<srcCS>& srcHess, 
    const PosDerivT<destCS, srcCS>& deriv, const PosDeriv2T<destCS, srcCS>& deriv2);

/** All-mighty routine for evaluating the value of a scalar function and its derivatives 
    in a different coordinate system (evalCS), and converting them to the target 
    coordinate system (outputCS). */
template<typename evalCS, typename outputCS>
EXP void evalAndConvert(const IScalarFunction<evalCS>& F,
    const PosT<outputCS>& pos,
    double* value=NULL,
    GradT<outputCS>* deriv=NULL,
    HessT<outputCS>* deriv2=NULL)
{
    bool needDeriv = deriv!=NULL || deriv2!=NULL;
    bool needDeriv2= deriv2!=NULL;
    GradT<evalCS> evalGrad;
    HessT<evalCS> evalHess;
    PosDerivT <outputCS, evalCS> coordDeriv;
    PosDeriv2T<outputCS, evalCS> coordDeriv2;
    const PosT<evalCS> evalPos = needDeriv ? 
        toPosDeriv<outputCS, evalCS>(pos, &coordDeriv, needDeriv2 ? &coordDeriv2 : 0) :
        toPos<outputCS, evalCS>(pos);
    // compute the function in transformed coordinates
    F.evalScalar(evalPos, value, needDeriv ? &evalGrad : 0, needDeriv2 ? &evalHess : 0);
    if(deriv)  // ... and convert gradient/hessian back to output coords if necessary.
        *deriv  = toGrad<evalCS, outputCS> (evalGrad, coordDeriv);
    if(deriv2)
        *deriv2 = toHess<evalCS, outputCS> (evalGrad, evalHess, coordDeriv, coordDeriv2);
}

/** The same routine for conversion of gradient and hessian 
    in the case that the computation requires the parameters of coordinate system evalCS */
template<typename evalCS, typename outputCS>
EXP void evalAndConvert(const IScalarFunction<evalCS>& F,
    const PosT<outputCS>& pos, const evalCS& coordsys,
    double* value=NULL, GradT<outputCS>* deriv=NULL, HessT<outputCS>* deriv2=NULL)
{
    bool needDeriv = deriv!=NULL || deriv2!=NULL;
    bool needDeriv2= deriv2!=NULL;
    GradT<evalCS> evalGrad;
    HessT<evalCS> evalHess;
    PosDerivT <outputCS, evalCS> coordDeriv;
    PosDeriv2T<outputCS, evalCS> coordDeriv2;
    const PosT<evalCS> evalPos = needDeriv ? 
        toPosDeriv<outputCS, evalCS>(pos, coordsys, &coordDeriv, needDeriv2 ? &coordDeriv2 : 0) :
        toPos<outputCS, evalCS>(pos, coordsys);
    // compute the function in transformed coordinates
    F.evalScalar(evalPos, value, needDeriv ? &evalGrad : 0, needDeriv2 ? &evalHess : 0);
    if(deriv)  // ... and convert gradient/hessian back to output coords if necessary.
        *deriv  = toGrad<evalCS, outputCS> (evalGrad, coordDeriv);
    if(deriv2)
        *deriv2 = toHess<evalCS, outputCS> (evalGrad, evalHess, coordDeriv, coordDeriv2);
}

/// trivial instantiation of the above function for the case that conversion is not necessary
template<typename CS> 
EXP void evalAndConvert(const IScalarFunction<CS>& F, const PosT<CS>& pos, 
    double* value, GradT<CS>* deriv, HessT<CS>* deriv2)
{  F.evalScalar(pos, value, deriv, deriv2); }

/** An even mightier routine for evaluating the value of a scalar function,
    its gradient and hessian, in a different coordinate system (evalCS), 
    and converting them to the target coordinate system (outputCS)
    through an intermediate coordinate system (intermedCS), 
    for the situation when a direct transformation is not available. */
template<typename evalCS, typename intermedCS, typename outputCS>
EXP void evalAndConvertTwoStep(const IScalarFunction<evalCS>& F,
    const PosT<outputCS>& pos,
    double* value=NULL,
    GradT<outputCS>* deriv=NULL,
    HessT<outputCS>* deriv2=NULL)
{
    bool needDeriv = deriv!=NULL || deriv2!=NULL;
    bool needDeriv2= deriv2!=NULL;
    GradT<evalCS> evalGrad;
    HessT<evalCS> evalHess;
    GradT<intermedCS> intermedGrad;
    HessT<intermedCS> intermedHess;
    PosDerivT <outputCS, intermedCS> coordDerivOI;
    PosDeriv2T<outputCS, intermedCS> coordDeriv2OI;
    PosDerivT <intermedCS, evalCS> coordDerivIE;
    PosDeriv2T<intermedCS, evalCS> coordDeriv2IE;
    const PosT<intermedCS> intermedPos = needDeriv ? 
        toPosDeriv<outputCS, intermedCS>(pos, &coordDerivOI, needDeriv2 ? &coordDeriv2OI : 0) :
        toPos<outputCS, intermedCS>(pos);
    const PosT<evalCS> evalPos = needDeriv ? 
        toPosDeriv<intermedCS, evalCS>(intermedPos, &coordDerivIE, needDeriv2 ? &coordDeriv2IE : 0) :
        toPos<intermedCS, evalCS>(intermedPos);
    // compute the function in transformed coordinates
    F.evalScalar(evalPos, value, needDeriv ? &evalGrad : 0, needDeriv2 ? &evalHess : 0);
    if(needDeriv)  // may be needed for either grad or hess (or both)
        intermedGrad = toGrad<evalCS, intermedCS> (evalGrad, coordDerivIE);
    if(deriv)
        *deriv  = toGrad<intermedCS, outputCS> (intermedGrad, coordDerivOI);
    if(deriv2) {
        intermedHess = toHess<evalCS, intermedCS> (evalGrad, evalHess, coordDerivIE, coordDeriv2IE);
        *deriv2 = toHess<intermedCS, outputCS> (intermedGrad, intermedHess, coordDerivOI, coordDeriv2OI);
    }
}

/** The same routine for the case that evalCS requires the parameters of coordinate system */
template<typename evalCS, typename intermedCS, typename outputCS>
EXP void evalAndConvertTwoStep(const IScalarFunction<evalCS>& F,
    const PosT<outputCS>& pos, const evalCS& coordsys,
    double* value=NULL, GradT<outputCS>* deriv=NULL, HessT<outputCS>* deriv2=NULL)
{
    bool needDeriv = deriv!=NULL || deriv2!=NULL;
    bool needDeriv2= deriv2!=NULL;
    GradT<evalCS> evalGrad;
    HessT<evalCS> evalHess;
    GradT<intermedCS> intermedGrad;
    HessT<intermedCS> intermedHess;
    PosDerivT <outputCS, intermedCS> coordDerivOI;
    PosDeriv2T<outputCS, intermedCS> coordDeriv2OI;
    PosDerivT <intermedCS, evalCS> coordDerivIE;
    PosDeriv2T<intermedCS, evalCS> coordDeriv2IE;
    const PosT<intermedCS> intermedPos = needDeriv ? 
        toPosDeriv<outputCS, intermedCS>(pos, &coordDerivOI, needDeriv2 ? &coordDeriv2OI : 0) :
        toPos<outputCS, intermedCS>(pos);
    const PosT<evalCS> evalPos = needDeriv ? 
        toPosDeriv<intermedCS, evalCS>(intermedPos, coordsys, &coordDerivIE, needDeriv2 ? &coordDeriv2IE : 0) :
        toPos<intermedCS, evalCS>(intermedPos, coordsys);
    // compute the function in transformed coordinates
    F.evalScalar(evalPos, value, needDeriv ? &evalGrad : 0, needDeriv2 ? &evalHess : 0);
    if(needDeriv)  // may be needed for either grad or hess (or both)
        intermedGrad = toGrad<evalCS, intermedCS> (evalGrad, coordDerivIE);
    if(deriv)
        *deriv  = toGrad<intermedCS, outputCS> (intermedGrad, coordDerivOI);
    if(deriv2) {
        intermedHess = toHess<evalCS, intermedCS> (evalGrad, evalHess, coordDerivIE, coordDeriv2IE);
        *deriv2 = toHess<intermedCS, outputCS> (intermedGrad, intermedHess, coordDerivOI, coordDeriv2OI);
    }
}

/** Specialized conversion routine for spherically-symmetric functions.
    Convert the derivatives of a simple function that only depends on the spherical radius 
    into gradients and hessians in a target coordinate system (outputCS). */
template<typename outputCS>
void evalAndConvertSph(const math::IFunction& F,
    const PosT<outputCS>& pos,
    double* value=NULL,
    GradT<outputCS>* deriv=NULL,
    HessT<outputCS>* deriv2=NULL);

///@}
/// \section Miscellaneous routines
///@{

/// compute the total angular momentum for a point in the given coordinate system CoordT
template<typename CoordT>
EXP double Ltotal(const PosVelT<CoordT> &p);

/// compute the z-component of angular momentum for a point in the given coordinate system CoordT
template<typename CoordT>
EXP double Lz(const PosVelT<CoordT> &p);

///@}
/// \section 3d rotations
///@{

/** construct a 3d rotation matrix specified by Euler angles.
    Let (x,y,z) be the source reference frame, and (X,Y,Z) be the rotated target frame.
    The first rotation by angle alpha about the z axis creates an intermediate reference frame
    (x',y',z'), where the axis x' points along the line of nodes of the overall transformation.
    The second rotation by angle beta about the x' axis tilts the (x',y') plane by angle beta,
    creating a second intermediate reference frame (x'', y'', z'').
    The third rotation by angle gamma about the z'' axis does not change the orientation of z'',
    hence the final axis Z is the same as z'', and beta is the angle between Z and z.
    The composition of three rotations is described by an orthogonal rotation matrix R,
    such that a point with coordinates (x,y,z) in the original reference frame will have
    coordinates (X,Y,Z) in the rotated frame, specified by
    \code
    | X |   | mat[0]  mat[1]  mat[2] |   | x |
    | Y | = | mat[3]  mat[4]  mat[5] | * | y |
    | Z |   | mat[6]  mat[7]  mat[8] |   | z |
    \endcode
    where mat is the flattened matrix R in a row-major order.
    In other words, the point remains geometrically fixed, only the reference frame changes
    and the coordinates of this point change accordingly (passive rotation).
    The inverse rotation is produced by a triplet of angles (-gamma,-beta,-alpha),
    and its rotation matrix is simply the transpose of the forward rotation matrix.
    \param[in]  alpha, beta, gamma are three Euler rotation angles.
    \param[out] mat  will contain 9 elements of a 3x3 rotation matrix in row-major order.
*/
EXP void makeRotationMatrix(double alpha, double beta, double gamma, double mat[9]);

/** transform a 3d vector in cartesian coordinates using the rotation matrix */
inline void transformVector(const double mat[9], const double vec[3], double result[3])
{
    result[0] = mat[0] * vec[0] + mat[1] * vec[1] + mat[2] * vec[2];
    result[1] = mat[3] * vec[0] + mat[4] * vec[1] + mat[5] * vec[2];
    result[2] = mat[6] * vec[0] + mat[7] * vec[1] + mat[8] * vec[2];
}



}  // namespace coord


