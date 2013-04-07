#include "mriCellMaterial.h"

// Constructor Implementation
MRICellMaterial::MRICellMaterial(double visc, double den){
  viscosity = visc;
  density = den;
}

MRICellMaterial::~MRICellMaterial()
{
}

