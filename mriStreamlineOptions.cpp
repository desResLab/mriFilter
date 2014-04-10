#include "mriStreamlineOptions.h"
#include "mriConstants.h"

// Constructor
MRIStreamlineOptions::MRIStreamlineOptions(double minX, double maxX, double minY, double maxY, double minZ, double maxZ)
{
  // Initialize Plane Type
  planeSlice = kPlaneXY;
  // Minimim and Maximum Coord of the Slices
  // Minimum Values
  minBoxCoords[0] = minX;
  minBoxCoords[1] = minY;
  minBoxCoords[2] = minZ;
  // Maximum Values
  maxBoxCoords[0] = maxX;
  maxBoxCoords[1] = maxY;
  maxBoxCoords[2] = maxZ;
  // Number Of Grid Points
  gridTotals[0] = 50;
  gridTotals[1] = 50;
  // Parameters
  deltaT = 0.1;
  totalT = 1000.0;
}

MRIStreamlineOptions::~MRIStreamlineOptions(){
}

// Set Limits
void MRIStreamlineOptions::setLimits(double minX, double maxX, double minY, double maxY, double minZ, double maxZ){
  // Minimum Values
  minBoxCoords[0] = minX;
  minBoxCoords[1] = minY;
  minBoxCoords[2] = minZ;
  // Maximum Values
  maxBoxCoords[0] = maxX;
  maxBoxCoords[1] = maxY;
  maxBoxCoords[2] = maxZ;
}



