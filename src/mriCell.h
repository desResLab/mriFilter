#ifndef MRICELL_H
#define MRICELL_H

# include "mriTypes.h"

class mriCell{
  public:
    double concentration;
    mriDoubleVec velocity;
    mriDoubleVec auxVector;
    // Constructor and Destructor
    mriCell();
    virtual ~mriCell();
    // Member Functions
    double getQuantity(int qtyID);
    void   setQuantity(int qtyID, double value);
};

#endif // MRICELL_H
