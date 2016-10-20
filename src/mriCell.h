#ifndef MRICELL_H
#define MRICELL_H

# include "mriTypes.h"

class MRICell{
  public:
    double concentration;
    MRIDoubleVec velocity;
    MRIDoubleVec auxVector;
    // Constructor and Destructor
    MRICell();
    ~MRICell();
    // Member Functions
    double getQuantity(int qtyID);
    void setQuantity(int qtyID, double value);
};

#endif // MRICELL_H
