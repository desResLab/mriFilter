#ifndef MRISCAN_H
#define MRISCAN_H
#include <stdio.h>
#include <string>
#include <vector>

#include "mriCell.h"
#include "mriVolData.h"
#include "mriCellMaterial.h"
#include "mriSamplingOptions.h"
#include "mriStreamline.h"
#include "mriStreamlineOptions.h"
#include "mriThresholdCriteria.h"
#include "mriConstants.h"
#include "mriImagedata.h"
#include "mriExpansion.h"
#include "mriOptions.h"

// =====================================================
// MAIN CLASS FOR BOTH UNSTRUCTURED AND STRUCTURED SCANS
// =====================================================
class MRIScan{
  public:
    // CONSTRUCTOR
    MRIScan();
    // ============
    // DATA MEMBERS
    // ============s
    // Domain Dimension
    double domainSizeMin[3];
    double domainSizeMax[3];
    // Envelope Velocities
    double maxVelModule;
    // Velocities And Concentrations for all Measure Points
    int totalCellPoints;
    std::vector<MRICell> cellPoints;
    // Utility Functions
    bool hasPressureGradient;
    bool hasRelativePressure;
    bool hasReynoldsStress;
    double scanTime;
    // MRI Expansion
    MRIExpansion* expansion;
    // ================
    // COMMON FUNCTIONS
    // ================
    // SCAN TRAVERSAL ALGORITHMS
    bool AreThereNotVisitedNeighbor(int cell, bool* visitedCell);
    bool AreThereVisitedNeighbor(int cell, bool* visitedCell, bool* isBoundaryCell, int &visitedNeighbor);
    int  findFirstNotVisite(int cellTotal, bool* visitedCell, std::vector<int> cellStack);

    // SAMPLING
    void SampleVelocities(MRISamplingOptions SamplingOptions);
    // PRESSURE    
    void EvalRelativePressure(int startingCell, double refPressure);
    void PerformPressureIterations();
    void EvalCellPressureGradients(int currentCell, MRICellMaterial material,
                                   double* timeDeriv, double** firstDerivs, double** secondDerivs,
                                   double** ReynoldsStressGrad,
                                   double* pressureGrad);
    int GetCellFromStack(std::vector<int> &cellStack, bool* visitedCell, bool* isBoundaryCell, bool &finished, bool& secondStage);
    // TURBULENCE MODELLING
    void EvalReynoldsStressComponent();
    void EvalPressureIterative(int currentCell, double currentValue, bool* visitedCell, std::vector<int> &cellStack,int& cellCount);
    // FILTERING
    void ApplyMedianFilter(int qtyID,int maxIt);

    // =================
    // VIRTUAL FUNCTIONS
    // =================
    // TOPOLOGY
    virtual void GetNeighbourCells(int CurrentCell, std::vector<int> &coords);
    virtual bool IsInnerCell(int Cell);
    // EVAL DERIVATIVES
    virtual void EvalSpaceDerivs(int currentCell, double** firstDerivs, double** secondDerivs){}
    // EVAL REYNOLDS STRESSES
    virtual void EvalReynoldsStressGradient(int currentCell, double** ReynoldsStressGradient){}    
    // COMPATIBILITY OF DIFFERENT SCANS
    virtual bool isCompatibleWith(MRIScan* secondScan);

};

#endif // MRISCAN_H
