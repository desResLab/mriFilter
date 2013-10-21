#ifndef MRIEXPANSION_H
#define MRIEXPANSION_H

#include <string>
#include <vector>

class MRIExpansion
{
public:
    // Data Members
    int totalVortices;
    double* constantFluxCoeff = nullptr;
    double* vortexCoeff = nullptr;
    // Constructor
    MRIExpansion(int totVortex);
    // Copy Constructor
    MRIExpansion(MRIExpansion* otherExp);
    // Other Constructor
    MRIExpansion(std::vector<double> Expansion);

    // Distructor
    ~MRIExpansion();
    // MEMBER FUNCTIONS
    // Threshold Vortex Expansion
    void ApplyVortexThreshold(double ratio);
    void WriteToFile(std::string outFile);
};

#endif // MRIEXPANSION_H
