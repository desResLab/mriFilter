#ifndef MRITOPOLOGY_H
#define MRITOPOLOGY_H

# include "mriTypes.h"
# include "mriUtils.h"
# include "mriIO.h"
# include "mriException.h"

// ================
// GENERIC TOPOLOGY
// ================
class mriTopology{
  public:
    // GENERIC TOPOLOGY
    // Domain Dimension
    mriDoubleVec domainSizeMin;
    mriDoubleVec domainSizeMax;
    // Velocities And Concentrations for all Measure Points
    int totalCells;
    // Cells Centroid Coordinates
    mriDoubleMat cellLocations;
    // Cells Topology
    mriIntMat cellConnections;
    mriIntMat cellFaces;
    // Face Topology
    mriIntMat faceCells;
    mriIntMat faceConnections;
    mriIntMat faceEdges;
    mriDoubleVec faceArea;
    mriDoubleMat faceNormal;
    // Edge Topology
    mriIntMat edgeConnections;
    mriIntMat edgeFaces;

    // STRUCTURED GRID TOPOLOGY
    // Cells Totals
    mriIntVec cellTotals;
    mriDoubleMat cellLengths;
    // Auxiliary
    mriDoubleMat auxNodesCoords;

    mriTopology();
    mriTopology(const mriIntVec& totals,
                const mriDoubleVec& lengthX,
                const mriDoubleVec& lengthY,
                const mriDoubleVec& lengthZ,
                const mriDoubleVec& minlimits,
                const mriDoubleVec& maxlimits);
    
    virtual ~mriTopology();

    // MEMBER FUNCTIONS

    // BUILDING A TOPOLOGY
    void readFromVTK_ASCII(string vtkFileName, vtkStructuredPointsOptionRecord& vtkOptions);
    void readFromPLT_ASCII(string pltFileName, pltOptionRecord& pltOptions);
    
    // CREATE FROM TEMPLATE
    void createFromTemplate(mriSamples sampleType, const mriDoubleVec& params);

    // CREATION FROM STRUCTURED GRIDS
    void createGridFromVTKStructuredPoints(const vtkStructuredPointsOptionRecord& opts);

    // MAPPING
    void   mapIndexToCoords(int index, mriIntVec& intCoords);
    int    mapCoordsToIndex(int i, int j, int k);
    int    getTotalAuxNodes();
    int    getTotalFaces();
    double getEdgeFaceVortexCoeff(int edgeID, int faceID);
    void   getEdgeDirection(int edgeID, mriDoubleVec& edgeDirVector);
    void   getEdgeToFaceDirection(int edgeID, int faceID, mriDoubleVec& edgeFaceVector);
    void   getFaceCenter(int faceID, mriDoubleVec& fc);
    void   getEdgeCenter(int edgeID, mriDoubleVec& ec);
    void   buildAuxNodesCoords();
    void   getAuxNodeCoordinates(int nodeNum, mriDoubleVec& pos);
    void   mapIndexToAuxNodeCoords(int index, mriIntVec& intCoords);
    void   mapAuxCoordsToPosition(const mriIntVec& auxCoords, mriDoubleVec& pos);
    void   buildCellConnections();
    void   buildFaceConnections();
    int    addToFaceConnections(const mriIntVec& faceIds, vector<vector<mriFace* > >& AuxFirstNodeFaceList);
    int    addToEdgeConnections(const mriIntVec& edgeIds,vector<vector<mriEdge*> >& AuxFirstNodeEdgeList);
    void   buildFaceCells();
    void   buildEdgeConnections();
    void   buildFaceAreasAndNormals();
    void   getExternalFaceNormal(int cellID, int localFaceID, mriDoubleVec& extNormal);
    void   mapCoordsToPosition(const mriIntVec& coords, bool addMeshMinima, mriDoubleVec& pos);
    int    getAdjacentFace(int globalNodeNumber /*Already Ordered Globally x-y-z*/, int AdjType);
    void   getNeighborVortexes(int cellNumber,int dim,mriIntVec& idx);

    // MANIPULATIONS TO TOPOLOGY
    void   scalePositions(const mriDoubleVec& origin, double factor);
    void   crop(const mriDoubleVec& limitBox, mriBoolVec& indexes);

    // CHECK COMPATIBILITY
    bool isCompatibleTopology(mriTopology* topo);
};

#endif // MRITOPOLOGY_H
