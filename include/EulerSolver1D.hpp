/*----------------------------------------------------------------------*\
            Solver for unsteady inviscid compressible flows
    
    Name of program
        ColverT-Euler1D

    Description
        Solves the time-dependent compressible Euler equations using the
        Finite Volume Method and the 3rd order TVD Runge-Kutta 
        time marching scheme.

    Riemann solvers
        Intercell fluxes can be computed using the HLLC, HLLE, AUSM or 
        Roe's Riemann solver. More Riemann solvers will be implemented 
        in the future.
    
    Flux limiter
        The minmod slope limiter function is to be implemented.

    References
        - Einfeldt, On Godunov-Type Methods for Gas Dynamics (1988)
        - Wendroff, A Two-Dimensional HLLE Riemann Solver and Associated
          Godunov-Type Difference Scheme for Gas Dynamics (1999)
        - Toro, The HLLC Riemann Solver (2019)
        - Toro, Riemann Solvers and Numerical Methods for Fluid
          Dynamics (2009)
\*----------------------------------------------------------------------*/

#pragma once

#include <iostream>
#include <fstream>
#include <sstream>

#include "InitialValueProblem.hpp"
#include "FunctionsEuler1D.hpp"

namespace ColverT {
    
    namespace OneDimension {

        class EulerSolver1D 
        {
            public:
                EulerSolver1D();
                EulerSolver1D(std::string meshFileName, std::string caseFileName);
                void RunCalculations();
            private:
                float m_LBCoordinate;
                float m_RBCoordinate;
                unsigned int m_NCells;
                float m_tMax;
                float m_CFL;
                float m_Gamma;
                unsigned int m_IterMax;
                unsigned int m_InitialValueProblem;
                unsigned int m_SpatialScheme;
                unsigned int m_RiemannSolver;
        };

    }

}