#ifndef MRICELLMATERIAL_H
#define MRICELLMATERIAL_H

class MRICellMaterial
{
public:
  double viscosity;
  double density;
  // Constructor and Destructor
  MRICellMaterial(double visc, double den);
  ~MRICellMaterial();
};

#endif // MRICELLMATERIAL_H
