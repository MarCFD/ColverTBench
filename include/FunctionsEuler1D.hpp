/*----------------------------------------------------------------------*\
                    Main functions of the 1D Euler solver

    Description
        Main functions used by the 1D Euler solver, including some
        Riemann solvers, the 3rd order Runge-Kutta scheme and spatial
        reconstruction methods.
\*----------------------------------------------------------------------*/

#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>

#include "InitialValueProblem.hpp"

namespace ColverT {

    namespace OneDimension {
    
        void ComputeCellCenters(float *&cellCenters, float LBCoordinate, float dx, unsigned int NCells);
        /** 
         * Implement initial conditions based on the selected Initial Value Problem
         * Case 1. Sod's shock tube problem
         * Case 2. Woordward-Colella interacting blastwaves problem
        **/        
        void ApplyInitialConditions(unsigned int initialValueProblem, float gamma, float L, float *&cellCenters, float *&CSV, unsigned int NRKStages, unsigned int NComponents, unsigned int NCellsTotal);

        inline double SoundSpeed(double gamma, double pressure, double density);
        /**
         * Compute primitive Euler variables PV from conservative Euler variables CSV
         * PV = (rho, u, P)
         * CSV = (rho, rho*u, E)
        */
        void ComputePVFromCSV(float *&PV, float *&CSV, unsigned int NCellsTotal, unsigned int NComponents, unsigned int NRKStages, float gamma);
        /**
         * Print to the console the main properties the user has selected for the simulation case
        **/    
        void PrintCase(unsigned int initialValueProblem, unsigned int spatialScheme, unsigned int RiemannSolver, unsigned int NCells, float tMax, float CFL);
        /**
         * Compute a stable time step by estimating the largest wave speed present throughout the domain at the current time level 
         * See Chapter 6 of the book written by Toro (2009)
        **/
        void ComputeStableTimeStep(float &dt, float gamma, float dx, float CFL, float t, float tMax, float*& CSV, unsigned int iter, unsigned int NComponents, unsigned int NCellsTotal, unsigned int NRKStages);

        void ApplyBoundaryConditions(float *&PV, float *&CSV, unsigned int RKStage, unsigned int NCellsTotal, unsigned int NComponents, unsigned int NRKStages);
        /**
         * Spatial reconstruction of the left and right states of the conservative Euler variables for each intercell
        **/
        void ReconstructLeftRightStates(float *&CSVLeftState, float *&CSVRightState, float *&CSV, unsigned int cell, unsigned int RKStage, unsigned int NComponents, unsigned int NRKStages);
        /**
         * HLLC Riemann solver as described in Toro (2009)
        **/
        void ComputeHLLCFluxes(float *&intercellFluxes, float *&localIntercellFlux, float *&CSVLeftState, float *&CSVRightState, float *&CSV, float gamma, unsigned int RKStage, unsigned int NIntercells, unsigned int NComponents, unsigned int NRKStages);
        /**
         * Switch function that calls the Riemann solver the user has chosen
        **/
        void ComputeIntercellFluxes(unsigned int RiemannSolver, float *&intercellFluxes, float *&localIntercellFlux, float *&CSVLeftState, float *&CSVRightState, float *&CSV, float gamma, unsigned int RKStage, unsigned int NIntercells, unsigned int NComponents, unsigned int NRKStages);
        /**
         * Update the flow fields to the next Runge-Kutta stage
         * In the last Runge-Kutta stage, the flow fields are attributed to the first Runge-Kutta stage for the next time iteration
        */
        void UpdateFlowFields(float *&PV, float *&CSV, float *&intercellFluxes, float gamma, float dt, float dx, unsigned int RKStage, unsigned int NCellsTotal, unsigned int NComponents, unsigned int NRKStages);

        void PerformRungeKutta3(float *&PV, float *&CSV, float *&intercellFluxes, float *&localIntercellFlux, float *&CSVLeftState, float *&CSVRightState, float gamma, float dt, float dx, unsigned int RiemannSolver, unsigned int NIntercells, unsigned int NCellsTotal, unsigned int NComponents, unsigned int NRKStages);
        /**
         * Write a results.dat file containing x-coordinates, density, velocity, pressure and temperature values at the output time tMax
        **/
        void WriteSolutionData(float *&PV, float *&cellCenters, unsigned int NCellsTotal, unsigned int NComponents, unsigned int NRKStages);

    }

}