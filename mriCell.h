#ifndef MIRCELL_H
#define MIRCELL_H

typedef double MRIReal;

class MRICell
{
	public:
    MRIReal concentration;
    MRIReal velocity[3];
    MRIReal position[3];
    MRIReal filteredVel[3];
    MRIReal pressGrad[3];
    MRIReal relPressure;
    // Constructor and Destructor
	  MRICell();
	  ~MRICell();
    // Member Functions
    double getQuantity(int qtyID);
};

#endif // MIRCELL_H
