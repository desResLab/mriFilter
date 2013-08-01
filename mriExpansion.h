#ifndef MRIEXPANSION_H
#define MRIEXPANSION_H

class MRIExpansion
{
public:
    // Data Members
    int totalVortices;
    double* constantFluxCoeff = nullptr;
    double* vortexCoeff = nullptr;
    // Constructor
    MRIExpansion(int totVortex);
    // Distructor
    ~MRIExpansion();
    // MEMBER FUNCTIONS
    // Threshold Vortex Expansion
    void ApplyVortexThreshold(double ratio);
};

#endif // MRIEXPANSION_H
