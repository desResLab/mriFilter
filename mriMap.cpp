#include <math.h>
#include <algorithm>

#include "mriScan.h"
#include "mriStructuredScan.h"
#include "mriConstants.h"
#include "mriException.h"
#include "mriUtils.h"

// =====================================
// CHECK IF IS BOUNDARY OR INTERNAL CELL
// =====================================
bool MRIScan::IsInnerCell(int Cell){
  // Init Result
  bool isInside = true; 
  std::vector<int> others;
	// Get Neighbors
  GetNeighbourCells(Cell,others);
  // If There Are Zero then False
  for(int loopA=0;loopA<k3DNeighbors;loopA++){
		isInside = ((isInside)&&(others[loopA]>-1));
	} 
  // RETURN
  return isInside;
}

// Get Unit Vector From Current Cell To Face Centre
void MRIScan::GetUnitVector(int CurrentCell, double* GlobalFaceCoords, double* &myVect){
  // Get Vector
  myVect[0] = GlobalFaceCoords[0]-cellPoints[CurrentCell].position[0];
  myVect[1] = GlobalFaceCoords[1]-cellPoints[CurrentCell].position[1];
  myVect[2] = GlobalFaceCoords[2]-cellPoints[CurrentCell].position[2];
  // Normalize
  MRIUtils::Normalize3DVector(myVect);
}

// ===================================
// GET GLOBAL COORDS FROM LOCAL COORDS
// ===================================
void MRIStructuredScan::GetGlobalCoords(int DimNumber, int SliceNumber, double FaceCoord1, double FaceCoord2, double* &globalCoords){
  // Sum up the slices
  double sliceValue = 0.0;
  for(int loopA=1;loopA<SliceNumber;loopA++){
    sliceValue += 0.5*(cellLengths[DimNumber][loopA-1] + cellLengths[DimNumber][loopA]);
  }
  switch(DimNumber){
    case 0:
      // X Direction
      globalCoords[0] = domainSizeMin[0] + sliceValue;
      globalCoords[1] = domainSizeMin[1] + FaceCoord1;
      globalCoords[2] = domainSizeMin[2] + FaceCoord2;
      break;
    case 1:
      // Y Direction
      globalCoords[0] = domainSizeMin[0] + FaceCoord1;
      globalCoords[1] = domainSizeMin[1] + sliceValue;
      globalCoords[2] = domainSizeMin[2] + FaceCoord2;
      break;
    case 2:
      // Z Direction
      globalCoords[0] = domainSizeMin[0] + FaceCoord1;
      globalCoords[1] = domainSizeMin[1] + FaceCoord2;
      globalCoords[2] = domainSizeMin[2] + sliceValue;
      break;
  }
}

/*// ===========================
// LOCAL TO GLOBAL FACE NUMBER
// ===========================
int MRIStructuredScan::FaceLocaltoGlobal(int LocalFace, int DimNumber, int SliceNumber){
  int cells1 = 0;
  double cellLength1 = 0.0;
  double cellLength2 = 0.0;
  // If Face is zero then return
  if(LocalFace == -1){
    return -1;
  }
  // Number Of Cells in Local 1
  switch(DimNumber){
    case 0:
      // YZ
      cells1 = cellTotals[1];
      cellLength1 = cellLength[1];
      cellLength2 = cellLength[2];
      break;
    case 1:
      // XZ
      cells1 = cellTotals[0];
      cellLength1 = cellLength[0];
      cellLength2 = cellLength[2];
      break;
    case 2:
      // XY
      cells1 = cellTotals[0];
      cellLength1 = cellLength[0];
      cellLength2 = cellLength[1];
      break;
  }
  // Find The Coords
  // CHECK IF IS OK!!!
  int iCoord = ((int)(LocalFace) / (int)(2*cells1+1));
  int jCoord = (     (LocalFace) %      (2*cells1+1));
  double faceCoord1,faceCoord2;
  if(jCoord>(cells1-1)){
    jCoord = jCoord - cells1;
    faceCoord1 = (jCoord)*cellLength1-0.5*cellLength1;
    faceCoord2 = (iCoord)*cellLength2;
  }else{
    faceCoord1 = (jCoord)*cellLength1;
    faceCoord2 = (iCoord)*cellLength2-0.5*cellLength2;
  }
  double* GlobalFaceCoords = new double[kNumberOfDimensions];
  // Pass to Global Coords
  GetGlobalCoords(DimNumber,SliceNumber,faceCoord1,faceCoord2,GlobalFaceCoords);
  // Find Global Cell
  int CurrentCell = GetCellNumber(GlobalFaceCoords);
  // Allocate
  double* UnitVector = new double[kNumberOfDimensions];
  // Find the Vector
  GetUnitVector(CurrentCell,GlobalFaceCoords,UnitVector);
  // Get Face Number
  int globalFace = GetFacewithCellVector(CurrentCell,UnitVector);
  // Deallocate
  delete [] UnitVector;
  delete [] GlobalFaceCoords;
  // Return Global Face Number
  return globalFace;
}*/

// Map To Cells Coords
void MRIStructuredScan::MapIndexToCoords(int index, int* intCoords){
  int CurrentIndex = index;
  intCoords[2] = (int)(CurrentIndex/(cellTotals[0]*cellTotals[1]));
  CurrentIndex = (CurrentIndex-intCoords[2]*cellTotals[0]*cellTotals[1]);
  intCoords[1] = (int)(CurrentIndex/cellTotals[0]);
  CurrentIndex = CurrentIndex-intCoords[1]*cellTotals[0];
  intCoords[0] = CurrentIndex;
}

// Map To Cells Coords
void MRIStructuredScan::MapIndexToAuxNodeCoords(int index, int* intCoords){
  int CurrentIndex = index;
  int totalX = cellTotals[0] + 1;
  int totalY = cellTotals[1] + 1;
  intCoords[2] = (int)(CurrentIndex/(totalX*totalY));
  CurrentIndex = (CurrentIndex-intCoords[2]*totalX*totalY);
  intCoords[1] = (int)(CurrentIndex/totalX);
  CurrentIndex = CurrentIndex-intCoords[1]*totalX;
  intCoords[0] = CurrentIndex;
}


// Map From Cells Coords
int MRIStructuredScan::MapCoordsToIndex(int i, int j, int k){
	// C++ INDEXES ZERO BASED
  return k*(cellTotals[0]*cellTotals[1])+j*(cellTotals[0])+i;
}

// Map Cell Number
int MRIStructuredScan::GetCellNumber(MRIReal* coords){
  // Check Indexes
  int i = MRIUtils::FindHowMany(coords[0]-domainSizeMin[0],cellLengths[0]);
  int j = MRIUtils::FindHowMany(coords[1]-domainSizeMin[1],cellLengths[1]);
  int k = MRIUtils::FindHowMany(coords[2]-domainSizeMin[2],cellLengths[2]);
  // If OutSide Project To Max Index
  if(i>(cellTotals[0]-1)){
    i = (cellTotals[0]-1);
  }
  if(j>(cellTotals[1]-1)){
    j = (cellTotals[1]-1);
  }
  if(k>(cellTotals[2]-1)){
    k = (cellTotals[2]-1);
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
  int index = MapCoordsToIndex(i,j,k);
  return index;
}

// Compute Global Permutation
void MRIStructuredScan::GetGlobalPermutation(std::vector<int> &GlobalPerm){
  MRIReal Norm = 0.0;
  MRIReal VelNorm = 0.0;
  int invalidCount = 0;
  // Fill Permutation Vector
  int NumberOfZeroNorms = 0;
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    Norm = MRIUtils::Do3DEucNorm(cellPoints[loopA].position);
    VelNorm = MRIUtils::Do3DEucNorm(cellPoints[loopA].velocity);
    if (Norm<kMathZero){
      invalidCount++;
    }
    if ((Norm<kMathZero)&&(VelNorm<kMathZero)&&(invalidCount>1)){
      GlobalPerm.push_back(-1);
      NumberOfZeroNorms++;
    }else{
      GlobalPerm.push_back(GetCellNumber(cellPoints[loopA].position));
    }
  }
}

// Get Local Adjacent Plane
int GetLocalAdjacentFace(int localNodeNumber, int totalX, int AdjType){
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
int MRIStructuredScan::GetAdjacentFace(int globalNodeNumber /*Already Ordered Globally x-y-z*/, int AdjType){
  // Get The Z Coord
  double currentZCoord = cellPoints[globalNodeNumber].position[2]-domainSizeMin[2];
  
	// Find The Node Number in The Current Plane
  int ZCompleteLevels = MRIUtils::FindHowMany(currentZCoord,cellLengths[2]);
  int localNodeNumber = globalNodeNumber-ZCompleteLevels*cellTotals[0]*cellTotals[1];

  // Find The Adjacent face in the Current Plane
  int localFaceNumber = GetLocalAdjacentFace(localNodeNumber,cellTotals[0],AdjType);

  // Map to Global Face Numbering: Differentiate if in Z
  if(AdjType == kfacePlusZ){
    return localFaceNumber+(ZCompleteLevels+1)*cellTotals[0]*cellTotals[1]+
                           (ZCompleteLevels+1)*(cellTotals[0]*(cellTotals[1]+1)+
                           (cellTotals[0]+1)*cellTotals[1]);
		
	}else if(AdjType == kfaceMinusZ){
    return localFaceNumber+(ZCompleteLevels)*cellTotals[0]*cellTotals[1]+
                           (ZCompleteLevels)*(cellTotals[0]*(cellTotals[1]+1)+
                           (cellTotals[0]+1)*cellTotals[1]);
		
	}else{
    return localFaceNumber+(ZCompleteLevels+1)*cellTotals[0]*cellTotals[1]+
                           (ZCompleteLevels)*(cellTotals[0]*(cellTotals[1]+1)+
                           (cellTotals[0]+1)*cellTotals[1]);
  }
}

// ================================
// GET VORTEXES ASSOCIATED TO CELLS
// ================================
void MRIStructuredScan::getNeighborVortexes(int cellNumber,int dim,int& idx1,int& idx2,int& idx3,int& idx4){
  int coords[3];
  int initialOffset = 0;
  MapIndexToCoords(cellNumber,coords);
  if(dim == 1){
    // YZ Plane
    initialOffset = coords[0]*((cellTotals[1]+1)*(cellTotals[2]+1)) + coords[2]*(cellTotals[1]+1);
    idx1 = initialOffset + coords[1];
    idx2 = idx1 + 1;
    idx3 = idx2 + (cellTotals[1]);
    idx4 = idx3 + 1;
  }else if(dim == 2){
    // XZ Plane
    initialOffset = (cellTotals[0]*(cellTotals[1]+1)*(cellTotals[2]+1)) + coords[1]*((cellTotals[0]+1)*(cellTotals[2]+1)) + coords[2]*(cellTotals[0]+1);
    idx1 = initialOffset + coords[0];
    idx2 = idx1 + 1;
    idx3 = idx2 + (cellTotals[0]);
    idx4 = idx3 + 1;
  }else{
    // XY Plane
    initialOffset = (cellTotals[0]*(cellTotals[1]+1)*(cellTotals[2]+1)) +
                    ((cellTotals[0]+1)*cellTotals[1]*(cellTotals[2]+1)) +
                    coords[2]*((cellTotals[0]+1)*(cellTotals[1]+1)) + coords[1]*(cellTotals[0]+1);
    idx1 = initialOffset + coords[0];
    idx2 = idx1 + 1;
    idx3 = idx2 + (cellTotals[0]);
    idx4 = idx3 + 1;
  }
}

// ==============================
// MAP INTEGER COORDS TO POSITION
// ==============================
void MRIStructuredScan::MapCoordsToPosition(int* coords, bool addMeshMinima, double* pos){
  // Loop on the three dimensions
  for(int loopA=0;loopA<3;loopA++){
    pos[loopA] = 0.0;
    for(int loopB=1;loopB<(coords[loopA]+1);loopB++){
      pos[loopA] += 0.5*(cellLengths[loopA][loopB-1] + cellLengths[loopA][loopB]);
    }
  }
  if(addMeshMinima){
    for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
      pos[loopA] += domainSizeMin[loopA];
    }
  }
}

// ========================================
// MAP AUXILIARY INTEGER COORDS TO POSITION
// ========================================
void MRIStructuredScan::MapAuxCoordsToPosition(int* auxCoords, double* pos){
  // Loop on the three dimensions
  for(int loopA=0;loopA<3;loopA++){
    pos[loopA] = 0.0;
    for(int loopB=0;loopB<auxCoords[loopA];loopB++){
      pos[loopA] += cellLengths[loopA][loopB];
    }
  }
  for(int loopA=0;loopA<3;loopA++){
    pos[loopA] = pos[loopA] + domainSizeMin[loopA] - 0.5*cellLengths[loopA][0];
  }
}


// ===================================
// CHECK IF STRUCTURED MESH IS UNIFORM
// ===================================
bool MRIStructuredScan::isUniform(){
  bool res = true;
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    // Assign Reference Length
    double refLength = cellLengths[0][0];
    // Loop through all the lengths to see of they are equal
    for(int loopB=0;loopB<cellLengths[loopA].size();loopB++){
      res = (res && (fabs(refLength-cellLengths[loopA][loopB])<kMathZero));
    }
  }
  return res;
}
