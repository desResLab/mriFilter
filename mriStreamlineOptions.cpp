#include "mriStreamlineOptions.h"
#include "mriConstants.h"

MRIStreamlineOptions::MRIStreamlineOptions()
{
}

MRIStreamlineOptions::~MRIStreamlineOptions()
{
}

void MRIStreamlineOptions::SetDefaultSLOptions(){
  planeSlice = kPlaneXY;
  // The Other Coord
  distanceFactor = 0.01;
  // Minimim and Maximum Coord of the Slices
  minCoordFactor[0] = 0.1;
  minCoordFactor[1] = 0.1;
  maxCoordFactor[0] = 0.9;
  maxCoordFactor[1] = 0.9;
  // Number Of Grid Points
  gridTotals[0] = 20;
  gridTotals[1] = 20;
  // Parameters
  deltaT = 0.1;
  totalT = 100.0;
}


