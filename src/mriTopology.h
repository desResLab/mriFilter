#ifndef MRITOPOLOGY_H
#define MRITOPOLOGY_H

# include "mriTypes.h"

// ================
// GENERIC TOPOLOGY
// ================
class MRITopology{
  public:
    // GENERIC TOPOLOGY
    // Domain Dimension
    double domainSizeMin[3];
    double domainSizeMax[3];
    // Velocities And Concentrations for all Measure Points
    int totalCells;
    // Cells Centroid Coordinates
    MRIDoubleMat cellLocations;
    // Cells Topology
    MRIIntMat cellConnections;
    MRIIntMat cellFaces;
    // Face Topology
    MRIIntMat faceCells;
    MRIIntMat faceConnections;
    MRIIntMat faceEdges;
    MRIDoubleVec faceArea;
    MRIDoubleMat faceNormal;
    // Edge Topology
    MRIIntMat edgeConnections;
    MRIIntMat edgeFaces;

    // STRUCTURED GRID TOPOLOGY
    // Cells Totals
    MRIIntVec cellTotals;
    MRIDoubleMat cellLengths;
    // Auxiliary
    MRIDoubleMat auxNodesCoords;

    MRITopology();
    virtual ~MRITopology();
};

#endif // MRITOPOLOGY_H
