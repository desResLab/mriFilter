#ifndef MRIEXPANSION_H
#define MRIEXPANSION_H

#include <string>
#include <vector>

class MRIExpansion
{
public:
    // Data Members
    int totalVortices;
    double* constantFluxCoeff;
    double* vortexCoeff;
    // Constructor
    MRIExpansion(int totVortex);
    // Copy Constructor
    MRIExpansion(MRIExpansion* otherExp);

    // Distructor
    ~MRIExpansion();
    // MEMBER FUNCTIONS
    // Fill From Vector
    void FillFromVector(std::vector<double> Expansion);
    // Threshold Vortex Expansion
    void ApplyVortexThreshold(int thresholdType, double ratio);
    // Eval 2-Norm of Coefficient Vector
    double Get2Norm(bool onlyVortex);
    // Write to text File
    void WriteToFile(std::string outFile);
};

#endif // MRIEXPANSION_H
