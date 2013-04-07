#ifndef MRISTREAMLINEOPTIONS_H
#define MRISTREAMLINEOPTIONS_H

class MRIStreamlineOptions
{
  public:
    // The Plane Of the Slice
    int planeSlice;
    // The Other Coord
    double distanceFactor;
    // Minimim and Maximum Coord of the Slices
    double minCoordFactor[2];
    double maxCoordFactor[2];
    // Number Of Grid Points
    int gridTotals[2];
    // Time Integration Parameters 
    double deltaT;
    double totalT;
    // Constructor and Destrustor
    MRIStreamlineOptions();
    ~MRIStreamlineOptions();
    
    // Member Functions
    void SetDefaultSLOptions();
};

#endif // MRISTREAMLINEOPTIONS_H
