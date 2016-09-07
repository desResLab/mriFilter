#ifndef MRISTREAMLINEOPTIONS_H
#define MRISTREAMLINEOPTIONS_H

#include "mriConstants.h"

class MRIStreamlineOptions
{
  public:
    // The Plane Of the Slice
    int planeSlice;
    // Use a Box to Define Streamlines
    double minBoxCoords[3];
    double maxBoxCoords[3];
    // Number Of Grid Points
    int gridTotals[2];
    // Time Integration Parameters 
    double deltaT;
    double totalT;
    // Constructor and Destrustor
    MRIStreamlineOptions(double minX, double maxX, double minY, double maxY, double minZ, double maxZ);
    ~MRIStreamlineOptions();
    // Member Functions
    void setLimits(double minX, double maxX, double minY, double maxY, double minZ, double maxZ);
    
};

#endif // MRISTREAMLINEOPTIONS_H
