#include <math.h>
#include "mriCell.h"
#include "mriConstants.h"
#include "mriException.h"

MRICell::MRICell(){
  // Initialize Cell
  concentration = 0.0;
  velocity.resize(3,0.0);
  auxVector.resize(3,0.0);
}

MRICell::~MRICell(){
}

// ============
// GET QUANTITY
// ============
double MRICell::getQuantity(int qtyID){
  switch(qtyID){
    case kQtyConcentration:
      return concentration;
      break;
    case kQtyVelocityX:
      return velocity[0];
      break;
    case kQtyVelocityY:
      return velocity[1];
      break;
    case kQtyVelocityZ:
      return velocity[2];
      break;
    default:
      throw MRIException("ERROR: Invalid quantity in MRICell::getQuantity.\n");
      break;
  }
}

// ============
// SET QUANTITY
// ============
void MRICell::setQuantity(int qtyID, double value){
  switch(qtyID){
    case kQtyConcentration:
      concentration = value;
      break;
    case kQtyVelocityX:
      velocity[0] = value;
      break;
    case kQtyVelocityY:
      velocity[1] = value;
      break;
    case kQtyVelocityZ:
      velocity[2] = value;
      break;
    default:
      throw MRIException("ERROR: Invalid quantity.\n");
      break;
  }
}



