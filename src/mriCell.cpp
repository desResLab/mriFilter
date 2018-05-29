#include <math.h>
#include "mriCell.h"
#include "mriConstants.h"
#include "mriException.h"

mriCell::mriCell(){
  // Initialize Cell
  concentration = 0.0;
  velocity.resize(3,0.0);
  auxVector.resize(3,0.0);
}

mriCell::~mriCell(){
}

// ============
// GET QUANTITY
// ============
double mriCell::getQuantity(int qtyID){
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
      throw mriException("ERROR: Invalid quantity in mriCell::getQuantity.\n");
      break;
  }
}

// ============
// SET QUANTITY
// ============
void mriCell::setQuantity(int qtyID, double value){
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
      throw mriException("ERROR: Invalid quantity.\n");
      break;
  }
}



