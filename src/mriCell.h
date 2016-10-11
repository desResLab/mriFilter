#ifndef MIRCELL_H
#define MIRCELL_H

class MRICell{
  public:
    double concentration;
    double velocity[3];
    double position[3];
    double auxVector[3];
    double ReStress[6];
    double pressGrad[3];
    double relPressure;
    // Constructor and Destructor
    MRICell();
    ~MRICell();
    // Member Functions
    double getQuantity(int qtyID);
    void setQuantity(int qtyID, double value);
};

#endif // MIRCELL_H
