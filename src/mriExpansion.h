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
    void fillFromVector(std::vector<double> Expansion);
    // Threshold Vortex Expansion
    void applyVortexThreshold(int thresholdType, double ratio);
    // Eval 2-Norm of Coefficient Vector
    double get2Norm(bool onlyVortex);
    // Write to text File
    void writeToFile(std::string outFile);
};

#endif // MRIEXPANSION_H
