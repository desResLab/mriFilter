#ifndef MRITURBMODEL_H
#define MRITURBMODEL_H

#include "mriScan.h"

// BASE CLASS FOR TURBULENCE MODELS
class MRITurbModel{ 
  public:
    // Default Constructor
    MRITurbModel();
    // Destructor
    virtual ~MRITurbModel();
    // VIRTUAL METHODS
    // Eval Turbulent Viscosity
    virtual void EvalTurbulentViscosity(MRIScan* currScan, double* nuT){}
    // Eval Turbulent Kinetic Energy
    virtual void EvalTurbulentKineticEnergy(MRIScan* currScan, double* kTurb){}
};

// Boussinesq Model
class MRITurbBoussinesq: public MRITurbModel{
  public:
    // Default Constructor
    MRITurbBoussinesq():MRITurbModel(){}
    // Distructor
    virtual ~MRITurbBoussinesq(){}
    // Eval Turbulent Viscosity    
    virtual void EvalTurbulentViscosity(MRIScan* currScan, double* nuT);
    // Eval Turbulent Kinetic Energy
    virtual void EvalTurbulentKineticEnergy(MRIScan* currScan, MRIThresholdCriteria* threshold, double* kTurb);
};

#endif // MRITURBMODEL_H
