#include "EulerSolver1D.hpp"

int main()
{
    ColverT::OneDimension::EulerSolver1D Solver("mesh1D.ini", "problem1D.ini");

    Solver.RunCalculations();

    return 0;
}