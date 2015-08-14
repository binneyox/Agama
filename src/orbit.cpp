#include "orbit.h"
#include <stdexcept>
#include <cmath>
#include "math_ode.h"

namespace orbit{

/** limit on the maximum number of steps in ODE solver */
static const unsigned int ODE_MAX_NUM_STEP = 1e6;


/* evaluate r.h.s. of ODE in different coordinate systems */
template<typename coordT>
class OrbitIntegrator: public math::IOdeSystem {
public:
    OrbitIntegrator(const potential::BasePotential& p) :
        potential(p) {};
    
    /** apply the equations of motion */
    virtual void eval(const double t, const math::OdeStateType& y, math::OdeStateType& dydt) const;
    
    /** return the size of ODE system - three coordinates and three velocities */
    virtual unsigned int size() const { return 6;}
    
    virtual bool isStdHamiltonian() const;
private:
    const potential::BasePotential& potential;
};

template<>
void OrbitIntegrator<coord::Car>::eval(const double /*t*/,
    const math::OdeStateType& y, math::OdeStateType& dydt) const
{
    coord::GradCar grad;
    potential.eval(coord::PosCar(y[0], y[1], y[2]), NULL, &grad);
    dydt[0] = y[3];
    dydt[1] = y[4];
    dydt[2] = y[5];  // time deriv. of position
    dydt[3] = -grad.dx;
    dydt[4] = -grad.dy;
    dydt[5] = -grad.dz;
}
template<>
bool OrbitIntegrator<coord::Car>::isStdHamiltonian() const { return true; }

template<>
void OrbitIntegrator<coord::Cyl>::eval(const double /*t*/,
    const math::OdeStateType& y, math::OdeStateType& dydt) const
{
    const coord::PosVelCyl p(&y.front());
    coord::GradCyl grad;
    potential.eval(p, NULL, &grad);
    double Rsafe = p.R!=0 ? p.R : 1e-100;  // avoid NAN in degenerate cases
    dydt[0] = p.vR;
    dydt[1] = p.vz;
    dydt[2] = p.vphi/Rsafe;
    dydt[3] = -grad.dR + pow_2(p.vphi) / Rsafe;
    dydt[4] = -grad.dz;
    dydt[5] = -(grad.dphi + p.vR*p.vphi) / Rsafe;
}
template<>
bool OrbitIntegrator<coord::Cyl>::isStdHamiltonian() const { return false; }
        
template<>
void OrbitIntegrator<coord::Sph>::eval(const double /*t*/,
    const math::OdeStateType& y, math::OdeStateType& dydt) const
{
    const coord::PosVelSph p(&y.front());
    coord::GradSph grad;
    potential.eval(p, NULL, &grad);
    double rsafe = p.r!=0 ? p.r : 1e-100;
    double sintheta = sin(p.theta);
    if(sintheta == 0) sintheta = 1e-100;
    double cottheta = cos(p.theta)/sintheta;
    dydt[0] = p.vr;
    dydt[1] = p.vtheta/rsafe;
    dydt[2] = p.vphi/(rsafe*sintheta);
    dydt[3] = -grad.dr + (pow_2(p.vtheta) + pow_2(p.vphi)) / rsafe;
    dydt[4] = (-grad.dtheta + pow_2(p.vphi)*cottheta - p.vr*p.vtheta) / rsafe;
    dydt[5] = (-grad.dphi/sintheta - (p.vr+p.vtheta*cottheta)*p.vphi) / rsafe;
}
template<>
bool OrbitIntegrator<coord::Sph>::isStdHamiltonian() const { return false; }

template<typename coordT>
class TrajectoryOutput {
public:
    TrajectoryOutput(const math::BaseOdeSolver& _solver,
        double _timeStep,
        std::vector<coord::PosVelT<coordT> >& _outputTraj) :
    solver(_solver), timeStep(_timeStep), outputTraj(_outputTraj), timePrev(0) {};
    void processStep() {
        double timeCurr = solver.getTime();
        int i1 = static_cast<int>(timePrev/timeStep);
        int i2 = static_cast<int>(timeCurr/timeStep);
        // store trajectory at regular intervals of time
        for(int i=i1; i<=i2; i++) {
            double t = timeStep*i;
            if(timePrev==0 || (t>timePrev && t<=timeCurr)) {
                double data[6];
                for(int c=0; c<6; c++)
                    data[c] = solver.value(t, c);
                outputTraj.push_back(coord::PosVelT<coordT>(data));
            }
        }
        timePrev = timeCurr;
    }
private:
    const math::BaseOdeSolver& solver;
    const double timeStep;
    std::vector<coord::PosVelT<coordT> >& outputTraj;
    double timePrev;
};

template<typename coordT>
unsigned int integrate(const potential::BasePotential& potential,
    const coord::PosVelT<coordT>& initialConditions,
    const double totalTime,
    const double outputTimestep,
    std::vector<coord::PosVelT<coordT> >& outputTrajectory,
    const double accuracy)
{
    if(outputTimestep <= 0) 
        throw std::invalid_argument("Orbit integration: output timestep is zero");
    OrbitIntegrator<coordT> odeSystem(potential);
    math::OdeStateType vars(odeSystem.size());
    initialConditions.unpack_to(&vars.front());
    math::BaseOdeSolver* solver = new math::OdeSolverDOP853(odeSystem, 0, accuracy);
    TrajectoryOutput<coordT> trajOutput(*solver, outputTimestep, outputTrajectory);
    try{
        solver->init(vars);
        bool finished = false;
        unsigned int numSteps = 0;
        while(!finished) {
            if(solver->step() <= 0 || numSteps >= ODE_MAX_NUM_STEP)  // signal of error
                finished = true;
            else {
                numSteps++;
                finished = solver->getTime() >= totalTime;
            }
            trajOutput.processStep();
        }
        delete solver;
        return numSteps;
    }
    catch(...) {
        delete solver;
        throw;
    }
}

// explicit template instantiations to make sure all of them get compiled
template unsigned int integrate(const potential::BasePotential&, const coord::PosVelCar&,
    const double, const double, std::vector<coord::PosVelCar>&, const double);
template unsigned int integrate(const potential::BasePotential&, const coord::PosVelCyl&,
    const double, const double, std::vector<coord::PosVelCyl>&, const double);
template unsigned int integrate(const potential::BasePotential&, const coord::PosVelSph&,
    const double, const double, std::vector<coord::PosVelSph>&, const double);

}  // namespace orbit