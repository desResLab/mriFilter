#ifndef MRITOPOLOGY_H
#define MRITOPOLOGY_H

# include "mriTypes.h"
# include "mriUtils.h"
# include "mriIO.h"
# include "mriException.h"

// ================
// GENERIC TOPOLOGY
// ================
class MRITopology{
  public:
    // GENERIC TOPOLOGY
    // Domain Dimension
    MRIDoubleVec domainSizeMin;
    MRIDoubleVec domainSizeMax;
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
    MRITopology(const MRIIntVec& totals,
                const MRIDoubleVec& lengthX,
                const MRIDoubleVec& lengthY,
                const MRIDoubleVec& lengthZ,
                const MRIDoubleVec& minlimits,
                const MRIDoubleVec& maxlimits);
    virtual ~MRITopology();

    // MEMBER FUNCTIONS

    // BUILDING A TOPOLOGY
    void readFromVTK_ASCII(string vtkFileName, vtkStructuredPointsOptionRecord& vtkOptions);
    void readFromPLT_ASCII(string pltFileName, pltOptionRecord& pltOptions);

    // CREATION FROM STRUCTURED GRIDS
    void createGridFromVTKStructuredPoints(const vtkStructuredPointsOptionRecord& opts);

    // MAPPING
    void   mapIndexToCoords(int index, MRIIntVec& intCoords);
    int    mapCoordsToIndex(int i, int j, int k);
    int    getTotalAuxNodes();
    int    getTotalFaces();
    double getEdgeFaceVortexCoeff(int edgeID, int faceID);
    void   getEdgeDirection(int edgeID, MRIDoubleVec& edgeDirVector);
    void   getEdgeToFaceDirection(int edgeID, int faceID, MRIDoubleVec& edgeFaceVector);
    void   getFaceCenter(int faceID, MRIDoubleVec& fc);
    void   getEdgeCenter(int edgeID, MRIDoubleVec& ec);
    void   buildAuxNodesCoords();
    void   getAuxNodeCoordinates(int nodeNum, MRIDoubleVec& pos);
    void   mapIndexToAuxNodeCoords(int index, MRIIntVec& intCoords);
    void   mapAuxCoordsToPosition(const MRIIntVec& auxCoords, MRIDoubleVec& pos);
    void   buildCellConnections();
    void   buildFaceConnections();
    int    addToFaceConnections(const MRIIntVec& faceIds, vector<vector<mriFace* > >& AuxFirstNodeFaceList);
    int    addToEdgeConnections(const MRIIntVec& edgeIds,vector<vector<mriEdge*> >& AuxFirstNodeEdgeList);
    void   buildFaceCells();
    void   buildEdgeConnections();
    void   buildFaceAreasAndNormals();
    void   getExternalFaceNormal(int cellID, int localFaceID, MRIDoubleVec& extNormal);
    void   mapCoordsToPosition(const MRIIntVec& coords, bool addMeshMinima, MRIDoubleVec& pos);
    int    getAdjacentFace(int globalNodeNumber /*Already Ordered Globally x-y-z*/, int AdjType);
    void   getNeighborVortexes(int cellNumber,int dim,MRIIntVec& idx);

    // MANIPULATIONS TO TOPOLOGY
    void   scalePositions(const MRIDoubleVec& origin, double factor);
    void   crop(const MRIDoubleVec& limitBox, MRIBoolVec& indexes);

    // CHECK COMPATIBILITY
    bool isCompatibleTopology(MRITopology* topo);
};

#endif // MRITOPOLOGY_H
