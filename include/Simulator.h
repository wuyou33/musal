/************
* Implements simulation loop for NMPC
* Created by Jean, April 8 2015
*******************************/

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <iostream>
#include <time.h>

#include "macros.h"
#include "ExplicitRKintegrator.h"
#include "NonlinOCPsolver.h"
#include "SpgOCPsolver.h"
#include "MuShoot.h"
#include "SolverOptions.h"
#include "MyTrajectory.h"

/**
** Implements simulation methods
**/
template <class dynT,class costT,class pconT,class mayerT> class Simulator {

public:
	/* Default constructor */
	Simulator();

	/* Destructor */	
	~Simulator();
		
	/* Copy constructor */
	Simulator(const Simulator<dynT,costT,pconT,mayerT>&);

	/* Assignement operator */
	Simulator<dynT,costT,pconT,mayerT> & operator= (const Simulator<dynT,costT,pconT,mayerT>&);
        
    /**
    ** Set methods
    **/
    /* Set kkt tolerance on parametric AL subproblem */
    inline void setTolKKTal(const double kktTolAL){
        this->kktTolAL = kktTolAL;
    }

    /* Set simulation time */
    inline void setSimulationTime(const double Tsim){
        assert(Tsim>POSZER);
        this->Tsim = Tsim;
    }

    /* Set simulation length */
    inline void setSimulationLength(const int lsim){
        assert(lsim>=1);
        this->lsim = lsim;
    }
    
    /* Set simulation sampling period */
    inline void setSamplingPeriod(const double ts){
        assert(ts>POSZER);
        this->ts = ts;
    }
    
    /* Set prediction time for NMPC controller */
    inline void setPredictionTime(const double Tpred){
        assert(Tpred>POSZER);
        mshoot.setPredictionTime(Tpred);
    }

    /* Set # shooting nodes for NMPC controller */
    inline void setNumShootingNodes(const int Ns){
        assert(Ns>=2);
        mshoot.setNumberShootingNodes(Ns);
    }
    
    /* Set simulation integrator step-size */
    inline void setStepSizeSimulation(const double hfix){
        assert(hfix>POSZER);
        s_hfix = hfix;
    }
    
    /* Set # integration steps in simulation integrator */
    inline void setNumIntStepSimulation(const int s_nsteps){
        assert(s_nsteps>=1);
        this->s_nsteps = s_nsteps;
    }

    /* Set multiple-shooting integrator step-size */
    inline void setStepSizeShooting(const double hfix){
        assert(hfix>POSZER);
        m_hfix = hfix;
    }

    /* Set # integration steps in multiple-shooting integrator */
    inline void setNumIntStepShooting(const int m_nsteps){
        assert(m_nsteps>=1);
        this->m_nsteps = m_nsteps;
    }

    /* Set dynamics */
    inline void setDynamics(const dynT& nlsys){
        this->nlsys = nlsys;
    }
    
    /* - Set bounds on state,
    - Compute scaling transformation if required  */
    inline void setStateBounds(double* xmin,double* xmax){
        mshoot.setXmin(xmin);
        mshoot.setXmax(xmax);
    }

    /* - Set bounds on input,
    - Compute scaling transformation if required */
    inline void setInputBounds(double* umin,double* umax){
        mshoot.setUmin(umin);
        mshoot.setUmax(umax);
    }

    /* Set bounds on terminal state */
    inline void setTerminalStateBounds(double* xminT,double* xmaxT){
        mshoot.setXminT(xminT);
        mshoot.setXmaxT(xmaxT);
    }

    /* Set initial state */
    void setInitialState(double*);
    
    /* Set integrator */
    inline void setSimulationIntegrator(const ExplicitRKintegrator<dynT>& simInteg){
        this->simInteg = simInteg;
    }
    
    /* Set shooter */
    inline void setShooter(const MuShoot<dynT,costT,pconT,mayerT>& shoot){
        nmpcSolver.setShooter(shoot);
//        spgSolver.setShooter(shoot);
    }
    
    /* Set Butcher tableau for simulation integrator */
    void setButcherSimulation(int,double*,double*,double*);
    
    /* Set Butcher tableau for shooting integrator */
    void setButcherShooting(int,double*,double*,double*);
    
    /* Set solver options */
    inline void setSolverOptions(const SolverOptions& options){
        this->options = options;
    }
    
    /**
    ** Processing methods 
    **/
    /* Initialisation method */
    void init();
    
    /* Simulation loop */
    void run();

#if defined(__unix__)
    /* To compute difference between 2 timespecs */
    timespec getTimeDiff(timespec,timespec);
#endif    

    /**
    * Read & save methods 
    */
    /* Read data from MAT files and store them in z & mu */
    void readMATfiles();

    /* Save closed-loop state trajectory */
    void saveStates(const char*);

    /* Save NMPC inputs */
    void saveInputs(const char*);

    /* Save simulation instants */
    void saveTimes(const char*);

    /* Save dimensions */
    void saveDimensions(const char*);

    /* Save KKT */
    void saveKKT(const char*);

    /* Save 2-norm of equality constraints */
    void saveFeas(const char*);

    /* Save solving time */
    void saveSolveTime(const char*);

    /* Save cumlative CG iterations */
    void saveCG(const char*);

    /**
    * For debugging
    */
#if SCAL_VAR
    /* Displays state transformation */
    void displayStateTransfo();

    /* Displays input transformation */
    void displayInputTransfo();
#endif

private:
    /**
    * User-defined objects
    */
    /* Nonlinear dynamics to be simulated in closed-loop with nmpcSolver */
    dynT nlsys;

    /* Reference trajectory */
    MyTrajectory trajRef;

    /**
    * Processing objects
    */
    /* Multiple-shooting object */
    MuShoot<dynT,costT,pconT,mayerT> mshoot;

    /* Trust-region augmented Lagrangian solver */
    NonlinOCPsolver<dynT,costT,pconT,mayerT> nmpcSolver;

    /* SPG solver */
    SpgOCPsolver<dynT,costT,pconT,mayerT> spgSolver;

    /* Simulation integrator */
    ExplicitRKintegrator<dynT> simInteg; 
    
    /* Solver options object */
    SolverOptions options;
    
    /**
     * Simulation parameters & data 
     */
    /* Allocation indicator,1 if xCL & uCL allocated, -1 otherwise */
    int allocated;
        
    /* KKT tolerance on parametric augmented Lagrangian problem */
    double kktTolAL;

    /* Simulation time */
    double Tsim;
    
    /* Simulation sampling period */
    double ts;

    /* Simulation length */
    int lsim;
    
    /* Closed-loop state trajectory, ALLOCATED */
    double* xCL;
    
    /* Closed-loop input, ALLOCATED */
    double* uCL;

    /* Simulation time, ALLOCATED */
    double* simTimes;

    /* Closed-loop KKT conditions on AL, ALLOCATED */
    double* kktCL;

    /* Closed-loop NMPC equality constraints, ALLOCATED */
    double* necCL;

    /* Closed-loop solving time, ALLOCATED */
    double* solveTimeCL;

    /* Closed-loop cumulative CG iterations, ALLOCATED */
    int* cgCL;

    /* Array of integration states, ALLOCATED */    
    double* integStates;
    
    /** 
     * Integrators parameters 
     */
    /* Parameters for simulation integrator, ALLOCATED */
    int s_s;
    double *s_A;
    double *s_b;
    double *s_c;
    double s_hfix;
    int s_nsteps;
    
    /* Parameters for shooting integrator, ALLOCATED */
    int m_s;
    double *m_A;
    double *m_b;
    double *m_c;
    double m_hfix;
    int m_nsteps; 

};
#include "Simulator.tpp"

#endif


