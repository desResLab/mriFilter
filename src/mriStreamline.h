#ifndef MRISTREAMLINE_H
#define MRISTREAMLINE_H

#include "mriConstants.h"

class MRIStreamline
{
  public:
    int totalPoints;
    std::vector<double> xCoords;  
    std::vector<double> yCoords;
    std::vector<double> zCoords;
    //  Constructor and Distructor
    MRIStreamline();
    ~MRIStreamline();
    // Member Functions
    bool EvalSLIntersection(MRIDirection dir, double fixedCoord, double* &intCoords);
    void EvalSLArrivalPointDistribution(int totalSL, std::vector<MRIStreamline> streamlines, MRIDirection dir, double minCoord, double maxCoord, int totalSlices, std::vector<double> sliceCenter, std::vector<double> sliceNormArrivals);
    // Printing to File
    void AppendToFile(int index,std::string fileName);
};

#endif // MRISTREAMLINE_H
