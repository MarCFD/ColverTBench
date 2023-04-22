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

#include "EulerSolver1D.hpp"

namespace ColverT {
    
    namespace OneDimension {

        EulerSolver1D::EulerSolver1D() 
        {
            m_LBCoordinate = 0.0;
            m_RBCoordinate = 1.0;
            m_NCells = 10;
            m_tMax = 1.0;
            m_CFL = 0.5;
            m_Gamma = 1.4;
            m_IterMax = 100;
            m_InitialValueProblem = 1;
            m_SpatialScheme = 1;
            m_RiemannSolver = 1;
        }

        EulerSolver1D::EulerSolver1D(std::string meshFileName, std::string caseFileName)
        {
            std::string line;

           // Read mesh file 
            std::ifstream inputMeshFile(meshFileName);
            if (!inputMeshFile.is_open()) 
            {
                std::cout << "ERROR: Mesh file could not be found.\n";
            }
            else 
            {
                std::cout << "Reading mesh file...\n";
                while (!inputMeshFile.eof()) 
                {
                    getline(inputMeshFile,line);
                    std::stringstream sst1(line);
                    sst1 >> m_LBCoordinate;

                    getline(inputMeshFile,line);
                    std::stringstream sst2(line);
                    sst2 >> m_RBCoordinate;

                    getline(inputMeshFile,line);
                    std::stringstream sst3(line);
                    sst3 >> m_NCells;
                }
                inputMeshFile.close();
                std::cout << "Mesh reading done.\n";
            }

           // Read case file
            std::ifstream inputCaseFile(caseFileName);
            if (!inputCaseFile.is_open())
            {
                std::cout << "ERROR: Case file could not be found.\n";
            }
            else
            {
                std::cout << "Reading case file...\n";
                while (!inputCaseFile.eof()) 
                {
                    getline(inputCaseFile,line);
                    std::stringstream sstTmax(line);
                    sstTmax >> m_tMax;

                    getline(inputCaseFile,line);
                    std::stringstream sstIterMax(line);
                    sstIterMax >> m_IterMax;
                
                    getline(inputCaseFile,line);
                    std::stringstream sstCFL(line);
                    sstCFL >> m_CFL;

                    getline(inputCaseFile,line);
                    std::stringstream sstGamma(line);
                    sstGamma >> m_Gamma;

                    getline(inputCaseFile,line);
                    std::stringstream sstIVP(line);
                    sstIVP >> m_InitialValueProblem;
                    
                    getline(inputCaseFile,line);
                    std::stringstream sstSpatialScheme(line);
                    sstSpatialScheme >> m_SpatialScheme;
                    
                    getline(inputCaseFile,line);
                    std::stringstream sstRiemannSolver(line);
                    sstRiemannSolver >> m_RiemannSolver;
                }
                inputCaseFile.close();
                std::cout << "Case reading done.\n";
            }
        }

        void EulerSolver1D::RunCalculations() 
        {
            std::cout << "Initialising the simulation...\n";       
            float L = m_RBCoordinate - m_LBCoordinate;
            float dx = L/m_NCells;
            float dt = 0.0f;
            unsigned int NGhostCells = 2;
            unsigned int NCellsTotal = m_NCells + NGhostCells;
            unsigned int NIntercells = m_NCells + 1;
            unsigned int NComponents = 3;
            unsigned int NRKStages = 3;

            float* CSV = new float[NCellsTotal*NComponents*NRKStages];
            float* PV = new float[NCellsTotal*NComponents*NRKStages];
            float* intercellFluxes = new float[NIntercells*NComponents*NRKStages];
            float* CSVLeftState = new float[NComponents];
            float* CSVRightState = new float[NComponents];
            float* localIntercellFlux = new float[NComponents];
            float* cellCenters = new float[m_NCells];

            ComputeCellCenters(cellCenters, m_LBCoordinate, dx, m_NCells);

            ApplyInitialConditions(m_InitialValueProblem, m_Gamma, L, cellCenters, CSV, NRKStages, NComponents, NCellsTotal);
            
            ComputePVFromCSV(PV, CSV, NCellsTotal, NComponents, NRKStages, m_Gamma);

           // Set flow time to zero
            float t = 0.0;
            std::cout << "Initialisation done.\n";

            PrintCase(m_InitialValueProblem, m_SpatialScheme, m_RiemannSolver, m_NCells, m_tMax, m_CFL);

            std::cout << "Starting calculations...\n";
           // Start time loop
            for (unsigned int iter = 1; iter <= m_IterMax; iter++) {
                ComputeStableTimeStep(dt, m_Gamma, dx, m_CFL, t, m_tMax, CSV, iter, NComponents, NCellsTotal, NRKStages);
                PerformRungeKutta3(PV, CSV, intercellFluxes, localIntercellFlux, CSVLeftState, CSVRightState, m_Gamma, dt, dx, m_RiemannSolver, NIntercells, NCellsTotal, NComponents, NRKStages);
               // Advance flow time
                t = t + dt;
               // Print to the console the current iteration and time step size
                std::cout << "Iter: " << iter << " Time level: " << t << "\n";
                std::cout << "Time step size: " << dt << "\n";
               // Verify if the max output time has been reached
                if (std::fabs(t - m_tMax)/m_tMax <= 1e-8) {
                    std::cout << "Output time reached!\n";
                    break;
                }
            }
            std::cout << "Calculations done.\n";

            WriteSolutionData(PV, cellCenters, NCellsTotal, NComponents, NRKStages);

           // Deallocate memory
            delete[] CSV;
            delete[] PV;
            delete[] intercellFluxes;
            delete[] CSVLeftState;
            delete[] CSVRightState;
            delete[] localIntercellFlux;
            delete[] cellCenters;
        }

    }

}