/*----------------------------------------------------------------------*\
                         Initial Value Problems

    Description
        Base class and children classes of Initial Value Problems.
        Each Initial Value Problem initialises the conservative Euler
        variables CSV to the corresponding values (initial conditions).
\*----------------------------------------------------------------------*/

#pragma once

namespace ColverT {

    namespace OneDimension {
        /**
         * Base class: Initial Value Problem
        **/
        class InitialValueProblem
        {
            public:
                InitialValueProblem();
                virtual void ApplyCellConditions(float gamma, float L, float*& cellCenters, float*& CSV, unsigned int NRKStages, unsigned int NComponents, unsigned int NCellsTotal);
            protected:
                float m_DensityLeft, m_DensityRight;
                float m_VelocityLeft, m_VelocityRight;
                float m_PressureLeft, m_PressureRight;
        };
        /**
         * Child class of InitialValueProblem: Sod's shock tube problem
        **/
        class SodShockTube : private InitialValueProblem
        {
            public:
                SodShockTube();
                void ApplyCellConditions(float gamma, float L, float*& cellCenters, float*& CSV, unsigned int NRKStages, unsigned int NComponents, unsigned int NCellsTotal) override;
        };
        /**
         * Child class of InitialValueProblem: Woodward-Colella interacting blastwaves problem
        */
        class WoodwardColellaBlastwaves : private InitialValueProblem
        {
            public:
                WoodwardColellaBlastwaves();
                void ApplyCellConditions(float gamma, float L, float*& cellCenters, float*& CSV, unsigned int NRKStages, unsigned int NComponents, unsigned int NCellsTotal) override;
            private:
                float m_DensityMiddle; 
                float m_PressureMiddle;
        };

    } 

}