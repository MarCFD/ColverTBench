/*----------------------------------------------------------------------*\
                         Initial Value Problems

    Description
        Main class and sub-classes of Initial Value Problems.
        Each Initial Value Problem initialises the conservative Euler
        variables CSV to the corresponding values (initial conditions).
\*----------------------------------------------------------------------*/

#include "InitialValueProblem.hpp"

namespace ColverT {

    namespace OneDimension {
        
        InitialValueProblem::InitialValueProblem()
        {
            m_DensityLeft = 0.0f;
            m_DensityRight = 0.0f;
            m_VelocityLeft = 0.0f;
            m_VelocityRight = 0.0f;
            m_PressureLeft = 0.0f;
            m_PressureRight = 0.0f;
        }

        void InitialValueProblem::ApplyCellConditions(float gamma, float L, float*& cellCenters, float*& CSV, unsigned int NRKStages, unsigned int NComponents, unsigned int NCells)
        {
            for (unsigned int RKStage = 0; RKStage < NRKStages; RKStage++) {
                for (unsigned int cell = 0; cell < NCells; cell++) {
                    if (cellCenters[cell] <= 0.5*L) {
                        CSV[(cell+1)*NComponents*NRKStages + 0*NRKStages + RKStage] = m_DensityLeft;
                        CSV[(cell+1)*NComponents*NRKStages + 1*NRKStages + RKStage] = m_DensityLeft*m_VelocityLeft;
                        CSV[(cell+1)*NComponents*NRKStages + 2*NRKStages + RKStage] = m_PressureLeft/(gamma-1) + 0.5*m_DensityLeft*m_VelocityLeft*m_VelocityLeft;
                    }
                    else {
                        CSV[(cell+1)*NComponents*NRKStages + 0*NRKStages + RKStage] = m_DensityRight;
                        CSV[(cell+1)*NComponents*NRKStages + 1*NRKStages + RKStage] = m_DensityRight*m_VelocityRight;
                        CSV[(cell+1)*NComponents*NRKStages + 2*NRKStages + RKStage] = m_PressureRight/(gamma-1) + 0.5*m_DensityRight*m_VelocityRight*m_VelocityRight;
                    }
                }
            }
        }

        SodShockTube::SodShockTube()
        {
            m_DensityLeft = 1.0f;
            m_VelocityLeft = 0.0f;
            m_PressureLeft = 1.0f;
            m_DensityRight = 0.125f;
            m_VelocityRight = 0.0f;
            m_PressureRight = 0.1f;
        }

        void SodShockTube::ApplyCellConditions(float gamma, float L, float*& cellCenters, float*& CSV, unsigned int NRKStages, unsigned int NComponents, unsigned int NCellsTotal)
        {
            for (unsigned int cell = 1; cell < NCellsTotal-1; cell++) {
                if (cellCenters[cell-1] <= 0.5f*L) {
                    CSV[cell*NComponents*NRKStages + 0*NRKStages + 0] = m_DensityLeft;
                    CSV[cell*NComponents*NRKStages + 1*NRKStages + 0] = m_DensityLeft*m_VelocityLeft;
                    CSV[cell*NComponents*NRKStages + 2*NRKStages + 0] = m_PressureLeft/(gamma-1.f) + 0.5f*m_DensityLeft*m_VelocityLeft*m_VelocityLeft;
                }
                else {
                    CSV[cell*NComponents*NRKStages + 0*NRKStages + 0] = m_DensityRight;
                    CSV[cell*NComponents*NRKStages + 1*NRKStages + 0] = m_DensityRight*m_VelocityRight;
                    CSV[cell*NComponents*NRKStages + 2*NRKStages + 0] = m_PressureRight/(gamma-1.f) + 0.5f*m_DensityRight*m_VelocityRight*m_VelocityRight;
                }
            }
        }

        WoodwardColellaBlastwaves::WoodwardColellaBlastwaves()
        {
            m_DensityLeft = 1.0f;
            m_VelocityLeft = 0.0f;
            m_PressureLeft = 1000.0f;
            m_DensityRight = 1.0f;
            m_VelocityRight = 0.0f;
            m_PressureRight = 100.0f;
            m_DensityMiddle = 1.0f;
            m_PressureMiddle = 0.001f;
        }

        void WoodwardColellaBlastwaves::ApplyCellConditions(float gamma, float L, float*& cellCenters, float*& CSV, unsigned int NRKStages, unsigned int NComponents, unsigned int NCellsTotal)
        {
            for (unsigned int cell = 1; cell < NCellsTotal-1; cell++) {
                if (cellCenters[cell] <= 0.1f*L) {
                    CSV[cell*NComponents*NRKStages + 0*NRKStages + 0] = m_DensityLeft;
                    CSV[cell*NComponents*NRKStages + 1*NRKStages + 0] = m_DensityLeft*m_VelocityLeft;
                    CSV[cell*NComponents*NRKStages + 2*NRKStages + 0] = m_PressureLeft/(gamma-1.f) + 0.5f*m_DensityLeft*m_VelocityLeft*m_VelocityLeft;
                }
                if (cellCenters[cell] > 0.1f*L || cellCenters[cell] < 0.8f*L) {
                    CSV[cell*NComponents*NRKStages + 0*NRKStages + 0] = m_DensityMiddle;
                    CSV[cell*NComponents*NRKStages + 1*NRKStages + 0] = 0.0f;
                    CSV[cell*NComponents*NRKStages + 2*NRKStages + 0] = m_PressureMiddle/(gamma-1.f);
                }
                if (cellCenters[cell] >= 0.8f*L) {
                    CSV[cell*NComponents*NRKStages + 0*NRKStages + 0] = m_DensityRight;
                    CSV[cell*NComponents*NRKStages + 1*NRKStages + 0] = m_DensityRight*m_VelocityRight;
                    CSV[cell*NComponents*NRKStages + 2*NRKStages + 0] = m_PressureRight/(gamma-1.f) + 0.5f*m_DensityRight*m_VelocityRight*m_VelocityRight;
                }
            }
        }

    }

}