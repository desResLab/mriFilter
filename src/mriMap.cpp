#include <math.h>
#include <algorithm>

#include "mriSequence.h"
#include "mriConstants.h"
#include "mriException.h"
#include "mriUtils.h"

using namespace std;

// =================================
// GET CARTESIAN NEIGHBORS OF A CELL
// =================================
void MRISequence::getCartesianNeighbourCells(int CurrentCell,MRIIntVec& cellNeighbors, bool addself){
  MRIIntVec coords(3);
  cellNeighbors.clear();
  if(addself){
    cellNeighbors.push_back(CurrentCell);
  }
  //Get The Coordinates of the Current Cell
  mapIndexToCoords(CurrentCell,coords);
  // Get Neighbor
  // coords[0]
  if ((coords[0]-1)>=0){
    cellNeighbors.push_back(mapCoordsToIndex(coords[0]-1,coords[1],coords[2]));
  }else{
    cellNeighbors.push_back(-1);
  }
  // coords[1]
  if((coords[0]+1)<topology->cellTotals[0]){
    cellNeighbors.push_back(mapCoordsToIndex(coords[0]+1,coords[1],coords[2]));
  }else{
    cellNeighbors.push_back(-1);
  }
  // coords[2]
  if((coords[1]-1)>=0){
    cellNeighbors.push_back(mapCoordsToIndex(coords[0],coords[1]-1,coords[2]));
  }else{
    cellNeighbors.push_back(-1);
  }
  // coords[3]
  if((coords[1]+1)<topology->cellTotals[1]){
    cellNeighbors.push_back(mapCoordsToIndex(coords[0],coords[1]+1,coords[2]));
  }else{
    cellNeighbors.push_back(-1);
  }
  // coords[4]
  if((coords[2]-1)>=0){
    cellNeighbors.push_back(mapCoordsToIndex(coords[0],coords[1],coords[2]-1));
  }else{
    cellNeighbors.push_back(-1);
  }
    // coords[5]
  if((coords[2]+1)<topology->cellTotals[2]){
    cellNeighbors.push_back(mapCoordsToIndex(coords[0],coords[1],coords[2]+1));
  }else{
    cellNeighbors.push_back(-1);
  }
  // SCRAMBLE VECTOR !!!
  //std::random_shuffle(&cellNeighbors[0], &cellNeighbors[5]);
}


// =====================================
// CHECK IF IS BOUNDARY OR INTERNAL CELL
// =====================================
bool MRISequence::isInnerCell(int cell){
  // Init Result
  bool isInside = true; 
  std::vector<int> others;
	// Get Neighbors
  getCartesianNeighbourCells(cell,others,false);
  // If There Are Zero then False
  for(int loopA=0;loopA<k3DNeighbors;loopA++){
    isInside = ((isInside)&&(others[loopA]>-1));
  }
  // RETURN
  return isInside;
}

// Get Unit Vector From Current Cell To Face Centre
void MRISequence::getUnitVector(int CurrentCell, const MRIDoubleVec& GlobalFaceCoords, MRIDoubleVec& myVect){
  // Get Vector
  myVect[0] = GlobalFaceCoords[0] - topology->cellLocations[CurrentCell][0];
  myVect[1] = GlobalFaceCoords[1] - topology->cellLocations[CurrentCell][1];
  myVect[2] = GlobalFaceCoords[2] - topology->cellLocations[CurrentCell][2];
  // Normalize
  MRIUtils::normalize3DVector(myVect);
}

// ===================================
// GET GLOBAL COORDS FROM LOCAL COORDS
// ===================================
void MRISequence::getGlobalCoords(int DimNumber, int SliceNumber, double FaceCoord1, double FaceCoord2, MRIDoubleVec& globalCoords){
  // Sum up the slices
  double sliceValue = 0.0;
  for(int loopA=1;loopA<SliceNumber;loopA++){
    sliceValue += 0.5*(topology->cellLengths[DimNumber][loopA-1] + topology->cellLengths[DimNumber][loopA]);
  }
  switch(DimNumber){
    case 0:
      // X Direction
      globalCoords[0] = topology->domainSizeMin[0] + sliceValue;
      globalCoords[1] = topology->domainSizeMin[1] + FaceCoord1;
      globalCoords[2] = topology->domainSizeMin[2] + FaceCoord2;
      break;
    case 1:
      // Y Direction
      globalCoords[0] = topology->domainSizeMin[0] + FaceCoord1;
      globalCoords[1] = topology->domainSizeMin[1] + sliceValue;
      globalCoords[2] = topology->domainSizeMin[2] + FaceCoord2;
      break;
    case 2:
      // Z Direction
      globalCoords[0] = topology->domainSizeMin[0] + FaceCoord1;
      globalCoords[1] = topology->domainSizeMin[1] + FaceCoord2;
      globalCoords[2] = topology->domainSizeMin[2] + sliceValue;
      break;
  }
}

// Map To Cells Coords
void MRISequence::mapIndexToCoords(int index, MRIIntVec& intCoords){
  int CurrentIndex = index;
  intCoords[2] = (int)(CurrentIndex/(topology->cellTotals[0] * topology->cellTotals[1]));
  CurrentIndex = (CurrentIndex-intCoords[2] * topology->cellTotals[0] * topology->cellTotals[1]);
  intCoords[1] = (int)(CurrentIndex / topology->cellTotals[0]);
  CurrentIndex = CurrentIndex-intCoords[1] * topology->cellTotals[0];
  intCoords[0] = CurrentIndex;
}

// Map To Aux Cells Coords
void MRISequence::mapIndexToAuxNodeCoords(int index, MRIIntVec& intCoords){
  int CurrentIndex = index;
  int totalX = topology->cellTotals[0] + 1;
  int totalY = topology->cellTotals[1] + 1;
  intCoords[2] = (int)(CurrentIndex/(totalX*totalY));
  CurrentIndex = (CurrentIndex-intCoords[2]*totalX*totalY);
  intCoords[1] = (int)(CurrentIndex/totalX);
  CurrentIndex = CurrentIndex-intCoords[1]*totalX;
  intCoords[0] = CurrentIndex;
}


// Map From Cells Coords
int MRISequence::mapCoordsToIndex(int i, int j, int k){
	// C++ INDEXES ZERO BASED
  return k*(topology->cellTotals[0]*topology->cellTotals[1])+j*(topology->cellTotals[0])+i;
}

// Map Cell Number
int MRISequence::getCellNumber(const MRIDoubleVec& coords){
  // Check Indexes
  int i = MRIUtils::findHowMany(coords[0] - topology->domainSizeMin[0],topology->cellLengths[0]);
  int j = MRIUtils::findHowMany(coords[1] - topology->domainSizeMin[1],topology->cellLengths[1]);
  int k = MRIUtils::findHowMany(coords[2] - topology->domainSizeMin[2],topology->cellLengths[2]);
  // If OutSide Project To Max Index
  if(i>(topology->cellTotals[0]-1)){
    i = (topology->cellTotals[0]-1);
  }
  if(j>(topology->cellTotals[1]-1)){
    j = (topology->cellTotals[1]-1);
  }
  if(k>(topology->cellTotals[2]-1)){
    k = (topology->cellTotals[2]-1);
  }
  // If OutSide Project to Min Index}
  if(i<0){
    i = 0;
  }
  if(j<0){
    j = 0;
  }
  if(k<0){
    k = 0;
  }
  // Return
  int index = mapCoordsToIndex(i,j,k);
  return index;
}

// Get Velocity Norm At Specific Cell for all scans: sort of Frobenius Norm
double MRISequence::getVelocityNormAtCell(int cell){
  double norm = 0.0;
  for(int loopA=0;loopA<sequence.size();loopA++){
    norm += sequence[loopA]->cells[loopA].velocity[0] * sequence[loopA]->cells[loopA].velocity[0];
    norm += sequence[loopA]->cells[loopA].velocity[1] * sequence[loopA]->cells[loopA].velocity[1];
    norm += sequence[loopA]->cells[loopA].velocity[2] * sequence[loopA]->cells[loopA].velocity[2];
  }
  return sqrt(norm);
}

// Compute Global Permutation
void MRISequence::getGlobalPermutation(MRIIntVec& GlobalPerm){
  double Norm = 0.0;
  double VelNorm = 0.0;
  int invalidCount = 0;
  // Fill Permutation Vector
  int NumberOfZeroNorms = 0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    Norm = MRIUtils::do3DEucNorm(topology->cellLocations[loopA]);
    // Get Velocity Norm Form One Cell In all Scans
    VelNorm = getVelocityNormAtCell(loopA);
    if (Norm<kMathZero){
      invalidCount++;
    }
    if ((Norm<kMathZero)&&(VelNorm<kMathZero)&&(invalidCount>1)){
      GlobalPerm.push_back(-1);
      NumberOfZeroNorms++;
    }else{
      GlobalPerm.push_back(getCellNumber(topology->cellLocations[loopA]));
    }
  }
}

// Get Local Adjacent Plane
int getLocalAdjacentFace(int localNodeNumber, int totalX, int AdjType){
  switch(AdjType){
    case kfacePlusX:  
		  return (((int)localNodeNumber)/((int)totalX)) * (2 * totalX + 1) + (localNodeNumber % totalX) + totalX + 1;
			break;
    case kfaceMinusX: 
		  return ((int)localNodeNumber/(int)totalX) * (2*totalX+1) + (localNodeNumber % totalX) + totalX;
    case kfacePlusY : 
		  return (((int)localNodeNumber/(int)totalX)+1) * (2*totalX+1) + (localNodeNumber % totalX);
    case kfaceMinusY: 
		  return ((int)localNodeNumber/(int)totalX) * (2*totalX+1) + (localNodeNumber % totalX);
    case kfacePlusZ:  
		  return localNodeNumber;
    case kfaceMinusZ: 
		  return localNodeNumber;
	}
  return -1;
}

// Get Global Adjacent Plane
int MRISequence::getAdjacentFace(int globalNodeNumber /*Already Ordered Globally x-y-z*/, int AdjType){
  // Get The Z Coord
  double currentZCoord = topology->cellLocations[globalNodeNumber][2] - topology->domainSizeMin[2];
  
	// Find The Node Number in The Current Plane
  int ZCompleteLevels = MRIUtils::findHowMany(currentZCoord,topology->cellLengths[2]);
  int localNodeNumber = globalNodeNumber - ZCompleteLevels * topology->cellTotals[0] * topology->cellTotals[1];

  // Find The Adjacent face in the Current Plane
  int localFaceNumber = getLocalAdjacentFace(localNodeNumber, topology->cellTotals[0], AdjType);

  // Map to Global Face Numbering: Differentiate if in Z
  if(AdjType == kfacePlusZ){
    return localFaceNumber+(ZCompleteLevels+1)*topology->cellTotals[0]*topology->cellTotals[1]+
                           (ZCompleteLevels+1)*(topology->cellTotals[0]*(topology->cellTotals[1]+1)+
                           (topology->cellTotals[0]+1)*topology->cellTotals[1]);
		
	}else if(AdjType == kfaceMinusZ){
    return localFaceNumber+(ZCompleteLevels)*topology->cellTotals[0]*topology->cellTotals[1]+
                           (ZCompleteLevels)*(topology->cellTotals[0]*(topology->cellTotals[1]+1)+
                           (topology->cellTotals[0]+1)*topology->cellTotals[1]);
		
	}else{
    return localFaceNumber+(ZCompleteLevels+1)*topology->cellTotals[0]*topology->cellTotals[1]+
                           (ZCompleteLevels)*(topology->cellTotals[0]*(topology->cellTotals[1]+1)+
                           (topology->cellTotals[0]+1)*topology->cellTotals[1]);
  }
}

// ================================
// GET VORTEXES ASSOCIATED TO CELLS
// ================================
void MRISequence::getNeighborVortexes(int cellNumber,int dim,MRIIntVec& idx){
  // Loop through the edges
  MRIIntVec ElEdgeList;
  double currEdgeDirVector[3];
  int currFace = 0;
  int currEdge = 0;
  for(int loopA=0;loopA<topology->cellFaces[cellNumber].size();loopA++){
    currFace = topology->cellFaces[cellNumber][loopA];
    for(int loopB=0;loopB<topology->faceEdges[currFace].size();loopB++){
      currEdge = topology->faceEdges[currFace][loopB];
      MRIUtils::insertInList(currEdge,ElEdgeList);
    }
  }
  // Find the Edges Aligned with the Selected Dimension
  idx.clear();
  int currDir = 0;
  for(int loopA=0;loopA<ElEdgeList.size();loopA++){
    getEdgeDirection(ElEdgeList[loopA],currEdgeDirVector);
    if((fabs(currEdgeDirVector[1])<kMathZero)&&(fabs(currEdgeDirVector[2])<kMathZero)){
      currDir = 0;
    }else if((fabs(currEdgeDirVector[0])<kMathZero)&&(fabs(currEdgeDirVector[2])<kMathZero)){
      currDir = 1;
    }else if((fabs(currEdgeDirVector[0])<kMathZero)&&(fabs(currEdgeDirVector[1])<kMathZero)){
      currDir = 2;
    }else{
      throw MRIException("ERROR: Invalid Edge Direction in getNeighborVortexes.\n");
    }
    if(currDir == dim){
      idx.push_back(ElEdgeList[loopA]);
    }
  }
}

// ==============================
// MAP INTEGER COORDS TO POSITION
// ==============================
void MRISequence::mapCoordsToPosition(const MRIIntVec& coords, bool addMeshMinima, MRIDoubleVec& pos){
  // Loop on the three dimensions
  for(int loopA=0;loopA<3;loopA++){
    pos[loopA] = 0.0;
    for(int loopB=1;loopB<(coords[loopA]+1);loopB++){
      pos[loopA] += 0.5*(topology->cellLengths[loopA][loopB-1] + topology->cellLengths[loopA][loopB]);
    }
  }
  if(addMeshMinima){
    for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
      pos[loopA] += topology->domainSizeMin[loopA];
    }
  }
}

// ========================================
// MAP AUXILIARY INTEGER COORDS TO POSITION
// ========================================
void MRISequence::mapAuxCoordsToPosition(const MRIIntVec& auxCoords, MRIDoubleVec& pos){
  // Loop on the three dimensions
  for(int loopA=0;loopA<3;loopA++){
    pos[loopA] = 0.0;
    for(int loopB=0;loopB<auxCoords[loopA];loopB++){
      pos[loopA] += topology->cellLengths[loopA][loopB];
    }
  }
  for(int loopA=0;loopA<3;loopA++){
    pos[loopA] = pos[loopA] + topology->domainSizeMin[loopA] - 0.5*topology->cellLengths[loopA][0];
  }
}

// ====================================
// GET COORDINATEDS FOR AUXILIARY NODES
// ====================================
void MRISequence::getAuxNodeCoordinates(int nodeNum, MRIDoubleVec& pos){
  MRIIntVec intAuxCoords(3,0.0);
  // Map To Integer Coordinates
  mapIndexToAuxNodeCoords(nodeNum,intAuxCoords);
  // Map To Spatial Position
  mapAuxCoordsToPosition(intAuxCoords,pos);
}

// ==================
// GET EDGE DIRECTION
// ==================
void MRISequence::getEdgeDirection(int edgeID, double* edgeDirVector){
  int node1 = 0;
  int node2 = 0;
  MRIDoubleVec node1Pos(3,0.0);
  MRIDoubleVec node2Pos(3,0.0);

  // Get The Two Nodes
  node1 = topology->edgeConnections[edgeID][0];
  node2 = topology->edgeConnections[edgeID][1];

  // Eval Auxiliary Node Coordinates
  node1Pos[0] = topology->auxNodesCoords[node1][0];
  node1Pos[1] = topology->auxNodesCoords[node1][1];
  node1Pos[2] = topology->auxNodesCoords[node1][2];
  node2Pos[0] = topology->auxNodesCoords[node2][0];
  node2Pos[1] = topology->auxNodesCoords[node2][1];
  node2Pos[2] = topology->auxNodesCoords[node2][2];

  // Get the versor
  edgeDirVector[0] = node1Pos[0] - node2Pos[0];
  edgeDirVector[1] = node1Pos[1] - node2Pos[1];
  edgeDirVector[2] = node1Pos[2] - node2Pos[2];
  double modulus = (edgeDirVector[0] * edgeDirVector[0] + 
                    edgeDirVector[1] * edgeDirVector[1] + 
                    edgeDirVector[2] * edgeDirVector[2]);

  edgeDirVector[0] = fabs(edgeDirVector[0]/modulus);
  edgeDirVector[1] = fabs(edgeDirVector[1]/modulus);
  edgeDirVector[2] = fabs(edgeDirVector[2]/modulus);
}
