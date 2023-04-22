/*----------------------------------------------------------------------*\
                    Functions for the 1D Euler solver

    Description
        Main functions used by the 1D Euler solver, including some
        Riemann solvers, the 3rd order Runge-Kutta scheme and spatial
        reconstruction methods.
\*----------------------------------------------------------------------*/

#include "FunctionsEuler1D.hpp"

namespace ColverT {

    namespace OneDimension {

        void ComputeCellCenters(float*& cellCenters, float LBCoordinate, float dx, unsigned int NCells)
        {
            cellCenters[0] = LBCoordinate + 0.5f*dx;
            
            for (unsigned int i = 0; i < NCells; i++) {
                cellCenters[i] = LBCoordinate + (i+1)*dx - 0.5f*dx;
            }
        }
        
        void ApplyInitialConditions(unsigned int initialValueProblem, float gamma, float L, float*& cellCenters, float*& CSV, unsigned int NRKStages, unsigned int NComponents, unsigned int NCellsTotal)
        {
            switch (initialValueProblem) {
                case 1: {
                    SodShockTube SodShockTubeProblem;
                    SodShockTubeProblem.ApplyCellConditions(gamma, L, cellCenters, CSV, NRKStages, NComponents, NCellsTotal);
                    break;
                }
                case 2: {
                    WoodwardColellaBlastwaves WCBlastwavesProblem;
                    WCBlastwavesProblem.ApplyCellConditions(gamma, L, cellCenters, CSV, NRKStages, NComponents, NCellsTotal);
                    break;
                }
                default: {
                    InitialValueProblem DefaultProblem;
                    DefaultProblem.ApplyCellConditions(gamma, L, cellCenters, CSV, NRKStages, NComponents, NCellsTotal);
                }
            }
        }

        inline double SoundSpeed(double gamma, double pressure, double density)
        {
            return std::sqrt(gamma*pressure/density);
        }

        void ComputePVFromCSV(float*& PV, float*& CSV, unsigned int NCellsTotal, unsigned int NComponents, unsigned int NRKStages, float gamma)
        {
            for (unsigned int cell = 1; cell < NCellsTotal-1; cell++) {
                PV[cell*NComponents*NRKStages + 0*NRKStages + 0] = CSV[cell*NComponents*NRKStages + 0*NRKStages + 0];
                PV[cell*NComponents*NRKStages + 1*NRKStages + 0] = CSV[cell*NComponents*NRKStages + 1*NRKStages + 0]/CSV[cell*NComponents*NRKStages + 0*NRKStages + 0];
                PV[cell*NComponents*NRKStages + 2*NRKStages + 0] = (gamma-1.0f)*(CSV[cell*NComponents*NRKStages + 2*NRKStages + 0] - 0.5f*CSV[cell*NComponents*NRKStages + 1*NRKStages + 0]*PV[cell*NComponents*NRKStages + 1*NRKStages + 0]);
            }
        }

        void PrintCase(unsigned int initialValueProblem, unsigned int spatialScheme, unsigned int RiemannSolver, unsigned int NCells, float tMax, float CFL)
        {
            std::cout << "===== Setup case =====\n";

            std::cout << "Number of cells: " << NCells << "\n";
            std::cout << "Maximum output time: " << tMax << "\n";
            std::cout << "Courant number: " << CFL << "\n";
            
            switch(initialValueProblem) {
                case 1: {
                    std::cout << "Initial Value Problem: Sod's shock tube problem\n";
                    break;
                }
                case 2: {
                    std::cout << "Initial Value Problem: Woodward-Colella interacting blastwaves problem\n";
                    break;
                }
                default: {
                    std::cout << "WARNING: Invalid Initial Value Problem! Default values have been applied.\n";
                }
            }

            switch (spatialScheme) {
                case 1: {
                    std::cout << "Spatial reconstruction scheme: 1st Order scheme\n";
                    break;
                }
                default: {
                    std::cout << "WARNING: Invalid spatial reconstruction scheme! 1st Order scheme has been applied.\n";
                }
            }

            switch (RiemannSolver) {
                case 1: {
                    std::cout << "Riemann solver: HLLC Riemann solver\n";
                    break;
                }
                default: {
                    std::cout << "WARNING: Invalid Riemann solver! HLLC Riemann solver is applied.\n";
                }
            }

            std::cout << "======================\n";
        }

        void ComputeStableTimeStep(float &dt, float gamma, float dx, float CFL, float t, float tMax, float*& CSV, unsigned int iter, unsigned int NComponents, unsigned int NCellsTotal, unsigned int NRKStages)
        {
                float soundSpeed = 0.0f;
                float maxWaveSpeed = 0.0f;
               // Estimate the max wave speed present throughout the domain at the current time level
                for (unsigned int cell = 0; cell < NCellsTotal; cell++) {
                    float velocity = CSV[cell*NComponents*NRKStages + 1*NRKStages + 0]/CSV[cell*NComponents*NRKStages + 0*NRKStages + 0];
                    float pressure = (gamma-1.0f)*(CSV[cell*NComponents*NRKStages + 2*NRKStages + 0] - 0.5f*CSV[cell*NComponents*NRKStages + 0*NRKStages + 0]*velocity);
                    soundSpeed = std::sqrt(gamma*pressure/CSV[cell*NComponents*NRKStages + 0*NRKStages + 0]);
                    maxWaveSpeed = std::max(maxWaveSpeed, soundSpeed + std::fabs(velocity));
                }
               // Reduce the time step size for the first 10 iterations to prevent instabilities
                if (iter <= 10) {
                    dt = std::min(0.1f*dx/maxWaveSpeed, tMax-t);
                }
                else {
                    dt = std::min(CFL*dx/maxWaveSpeed, tMax-t);
                }
        }

        void ApplyBoundaryConditions(float *&PV, float *&CSV, unsigned int RKStage, unsigned int NCellsTotal, unsigned int NComponents, unsigned int NRKStages)
        {
           // Apply transmissive boundary conditions
           // TODO: implement reflective boundary conditions
            for (unsigned int component = 0; component < NComponents; component++) {
               // Left boundary
                CSV[0*NComponents*NRKStages + component*NRKStages + RKStage] = CSV[1*NComponents*NRKStages + component*NRKStages + RKStage];
                PV[0*NComponents*NRKStages + component*NRKStages + RKStage] = PV[1*NComponents*NRKStages + component*NRKStages + RKStage];
               // Right Boundary
                CSV[(NCellsTotal-1)*NComponents*NRKStages + component*NRKStages + RKStage] = CSV[(NCellsTotal-2)*NComponents*NRKStages + component*NRKStages + RKStage];
                PV[(NCellsTotal-1)*NComponents*NRKStages + component*NRKStages + RKStage] = PV[(NCellsTotal-2)*NComponents*NRKStages + component*NRKStages + RKStage];
            }
        }

        void ReconstructLeftRightStates(float *&CSVLeftState, float *&CSVRightState, float *&CSV, unsigned int cell, unsigned int RKStage, unsigned int NComponents, unsigned int NRKStages)
        {
           // TODO: Higher-order TVD reconstruction schemes
            for (unsigned int component = 0; component < NComponents; component++) {
                CSVLeftState[component] = CSV[cell*NComponents*NRKStages + component*NRKStages + RKStage];
                CSVRightState[component] = CSV[(cell+1)*NComponents*NRKStages + component*NRKStages + RKStage];
            }
        }

        void ComputeHLLCFluxes(float *&intercellFluxes, float *&localIntercellFlux, float *&CSVLeftState, float *&CSVRightState, float *&CSV, float gamma, unsigned int RKStage, unsigned int NIntercells, unsigned int NComponents, unsigned int NRKStages)
        {
            for (unsigned int cell = 0; cell < NIntercells; cell++) {
                ReconstructLeftRightStates(CSVLeftState, CSVRightState, CSV, cell, RKStage, NComponents, NRKStages);
               // Step 1: pressure estimate
               // PVRS scheme
                float uL = CSVLeftState[1]/CSVLeftState[0];
                float uR = CSVRightState[1]/CSVRightState[0];
                float pL = (gamma-1)*(CSVLeftState[2] - 0.5f*CSVLeftState[1]*uL);
                float pR = (gamma-1)*(CSVRightState[2] - 0.5f*CSVRightState[1]*uR);
                float aL = std::sqrt(gamma*pL/CSVLeftState[0]);
                float aR = std::sqrt(gamma*pR/CSVRightState[0]);
                float rhoBar = 0.5f*(CSVLeftState[0] + CSVRightState[0]);
                float aBar = 0.5f*(aL + aR);
                float pPVRS = 0.5f*(pL + pR) - 0.5f*(uR - uL)*rhoBar*aBar;
                
                float pStar = std::max(0.0f, pPVRS);
               // Step 2: wave speed estimates
               // Pressure-Based Wave Speed Estimate
                float qL = 0.0f;
                float qR = 0.0f;
                if (pStar <= pL) {
                    qL = 1.0f;
                }    
                else {
                    qL = std::sqrt(1 + (gamma+1)/(2*gamma)*((pStar/pL) - 1));
                }

                if (pStar <= pR) {
                    qR = 1.0f;
                }
                else {
                    qR = std::sqrt(1 + (gamma+1)/(2*gamma)*((pStar/pR) - 1));
                }
                float SL = uL - aL*qL;
                float SR = uR + aR*qR;
                float SStar = (pR - pL + CSVLeftState[0]*uL*(SL - uL) - CSVRightState[0]*uR*(SR - uR))/(CSVLeftState[0]*(SL - uL) - CSVRightState[0]*(SR - uR));           
               // Step 3: Compute HLLC flux
                float FL[3];
                float FR[3];
                float FSL[3];
                float FSR[3];
                float USL[3];
                float USR[3];
                FL[0] = CSVLeftState[1];
                FL[1] = CSVLeftState[1]*uL + pL;
                FL[2] = (CSVLeftState[2] + pL)*uL;
                FR[0] = CSVRightState[1];
                FR[1] = CSVRightState[1]*uR + pR;
                FR[2] = (CSVRightState[2] + pR)*uR;
                USL[0] = CSVLeftState[0]*(SL - uL)/(SL - SStar);
                USL[1] = CSVLeftState[0]*(SL - uL)/(SL - SStar)*SStar;
                USL[2] = CSVLeftState[0]*(SL - uL)/(SL - SStar)*(CSVLeftState[2]/CSVLeftState[0] + (SStar - uL)*(SStar + pL/(CSVLeftState[0]*(SL - uL))));
                USR[0] = CSVRightState[0]*(SR - uR)/(SR - SStar);
                USR[1] = CSVRightState[0]*(SR - uR)/(SR - SStar)*SStar;
                USR[2] = CSVRightState[0]*(SR - uR)/(SR - SStar)*(CSVRightState[2]/CSVRightState[0] + (SStar - uR)*(SStar + pR/(CSVRightState[0]*(SR - uR))));
                
                for (unsigned int component = 0; component < NComponents; component++) {
                    FSL[component] = FL[component] + SL*(USL[component] - CSVLeftState[component]);
                    FSR[component] = FR[component] + SR*(USR[component] - CSVRightState[component]);
                }

                if (SL >= 0) {
                    for (unsigned int component = 0; component < NComponents; component++) {
                        localIntercellFlux[component] = FL[component];
                    }
                }
                else if (SL <= 0 && SStar >= 0) {
                    for (unsigned int component = 0; component < NComponents; component++) {
                        localIntercellFlux[component] = FSL[component];
                    }
                }
                else if (SR >= 0 && SStar <= 0) {
                    for (unsigned int component = 0; component < NComponents; component++) {
                        localIntercellFlux[component] = FSR[component];
                    }
                }
                else if (SR <= 0) {
                    for (unsigned int component = 0; component < NComponents; component++) {
                        localIntercellFlux[component] = FR[component];
                    }
                }

                for (unsigned int component = 0; component < NComponents; component++) {
                    intercellFluxes[cell*NComponents*NRKStages + component*NRKStages + RKStage] = localIntercellFlux[component];
                }
            }
        }

        void ComputeIntercellFluxes(unsigned int RiemannSolver, float *&intercellFluxes, float *&localIntercellFlux, float *&CSVLeftState, float *&CSVRightState, float *&CSV, float gamma, unsigned int RKStage, unsigned int NIntercells, unsigned int NComponents, unsigned int NRKStages)
        {
            switch (RiemannSolver) {
                case 1: {
                    ComputeHLLCFluxes(intercellFluxes, localIntercellFlux, CSVLeftState, CSVRightState, CSV, gamma, RKStage, NIntercells, NComponents, NRKStages);
                    break;
                }
                default: {
                    ComputeHLLCFluxes(intercellFluxes, localIntercellFlux, CSVLeftState, CSVRightState, CSV, gamma, RKStage, NIntercells, NComponents, NRKStages);
                }
            }
        }

        void UpdateFlowFields(float *&PV, float *&CSV, float *&intercellFluxes, float gamma, float dt, float dx, unsigned int RKStage, unsigned int NCellsTotal, unsigned int NComponents, unsigned int NRKStages)
        {
            for (unsigned int cell = 1; cell < NCellsTotal-1; cell++) {
                switch(RKStage) {
                    case 0: {
                        for (unsigned int component = 0; component < NComponents; component++) {
                            CSV[cell*NComponents*NRKStages + component*NRKStages + 1] = CSV[cell*NComponents*NRKStages + component*NRKStages + 0] - (dt/dx)*(intercellFluxes[cell*NComponents*NRKStages + component*NRKStages + 0] - intercellFluxes[(cell-1)*NComponents*NRKStages + component*NRKStages + 0]);
                        }

                        PV[cell*NComponents*NRKStages + 0*NRKStages + 1] = CSV[cell*NComponents*NRKStages + 0*NRKStages + 1];
                        PV[cell*NComponents*NRKStages + 1*NRKStages + 1] = CSV[cell*NComponents*NRKStages + 1*NRKStages + 1]/CSV[cell*NComponents*NRKStages + 0*NRKStages + 1];
                        PV[cell*NComponents*NRKStages + 2*NRKStages + 1] = (gamma-1)*(CSV[cell*NComponents*NRKStages + 2*NRKStages + 1] - 0.5f*CSV[cell*NComponents*NRKStages + 1*NRKStages + 1]*PV[cell*NComponents*NRKStages + 1*NRKStages + 1]);

                        break;
                    }
                    case 1: {
                        for (unsigned int component = 0; component < NComponents; component++) {
                            CSV[cell*NComponents*NRKStages + component*NRKStages + 2] = 0.75f*CSV[cell*NComponents*NRKStages + component*NRKStages + 0] + 0.25f*CSV[cell*NComponents*NRKStages + component*NRKStages + 1] - (0.25f*dt/dx)*(intercellFluxes[cell*NComponents*NRKStages + component*NRKStages + 1] - intercellFluxes[(cell-1)*NComponents*NRKStages + component*NRKStages + 1]);
                        }

                        PV[cell*NComponents*NRKStages + 0*NRKStages + 2] = CSV[cell*NComponents*NRKStages + 0*NRKStages + 2];
                        PV[cell*NComponents*NRKStages + 1*NRKStages + 2] = CSV[cell*NComponents*NRKStages + 1*NRKStages + 2]/CSV[cell*NComponents*NRKStages + 0*NRKStages + 2];
                        PV[cell*NComponents*NRKStages + 2*NRKStages + 2] = (gamma-1)*(CSV[cell*NComponents*NRKStages + 2*NRKStages + 2] - 0.5f*CSV[cell*NComponents*NRKStages + 1*NRKStages + 2]*PV[cell*NComponents*NRKStages + 1*NRKStages + 2]);

                        break;
                    }
                    case 2: {
                        for (unsigned int component = 0; component < NComponents; component++) {
                            CSV[cell*NComponents*NRKStages + component*NRKStages + 0] = (1.f/3.f)*CSV[cell*NComponents*NRKStages + component*NRKStages + 0] + (2.f/3.F)*CSV[cell*NComponents*NRKStages + component*NRKStages + 2] - (2.f/3.f*dt/dx)*(intercellFluxes[cell*NComponents*NRKStages + component*NRKStages + 2] - intercellFluxes[(cell-1)*NComponents*NRKStages + component*NRKStages + 2]);
                        }

                        PV[cell*NComponents*NRKStages + 0*NRKStages + 0] = CSV[cell*NComponents*NRKStages + 0*NRKStages + 0];
                        PV[cell*NComponents*NRKStages + 1*NRKStages + 0] = CSV[cell*NComponents*NRKStages + 1*NRKStages + 0]/CSV[cell*NComponents*NRKStages + 0*NRKStages + 0];
                        PV[cell*NComponents*NRKStages + 2*NRKStages + 0] = (gamma-1)*(CSV[cell*NComponents*NRKStages + 2*NRKStages + 0] - 0.5f*CSV[cell*NComponents*NRKStages + 1*NRKStages + 0]*PV[cell*NComponents*NRKStages + 1*NRKStages + 0]);

                        break;
                    }
                }
            }
        }

        void PerformRungeKutta3(float *&PV, float *&CSV, float *&intercellFluxes, float *&localIntercellFlux, float *&CSVLeftState, float *&CSVRightState, float gamma, float dt, float dx, unsigned int RiemannSolver, unsigned int NIntercells, unsigned int NCellsTotal, unsigned int NComponents, unsigned int NRKStages)
        {
            for (unsigned int RKStage = 0; RKStage < NRKStages; RKStage++) {
                ApplyBoundaryConditions(PV, CSV, RKStage, NCellsTotal, NComponents, NRKStages);

                ComputeIntercellFluxes(RiemannSolver, intercellFluxes, localIntercellFlux, CSVLeftState, CSVRightState, CSV, gamma, RKStage, NIntercells, NComponents, NRKStages);

                UpdateFlowFields(PV, CSV, intercellFluxes, gamma, dt, dx, RKStage, NCellsTotal, NComponents, NRKStages);
            }
        }

        void WriteSolutionData(float *&PV, float *&cellCenters, unsigned int NCellsTotal, unsigned int NComponents, unsigned int NRKStages)
        {
            std::cout << "Writing solutions...\n";
            FILE *dataFile;
            dataFile = fopen("results.dat", "w");

            if (dataFile == NULL) {
                printf("\nERROR when opening file!\n");
            }

            fprintf(dataFile, "VARIABLES=\"X\",\"RHO\",\"U\",\"P\",\"T\",\n");

            for (unsigned int cell = 1; cell < NCellsTotal-1; cell++) 
            {
                fprintf(dataFile, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", cellCenters[cell-1], PV[cell*NComponents*NRKStages + 0*NRKStages + 0], PV[cell*NComponents*NRKStages + 1*NRKStages + 0], PV[cell*NComponents*NRKStages + 2*NRKStages + 0], PV[cell*NComponents*NRKStages + 2*NRKStages + 0]/PV[cell*NComponents*NRKStages + 0*NRKStages + 0]);
            }
            
            fclose(dataFile);

            std::cout << "Solutions are written in results.dat!\n";
        }

    }

}
