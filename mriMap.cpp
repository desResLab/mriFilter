#include <math.h>
#include "mriScan.h"
#include "mriConstants.h"
#include "mriException.h"
#include "mriUtils.h"

// Check If Cell Is Inside or On the Boundary
bool MRIScan::IsInnerCell(int Cell){
  // Init Result
  bool isInside = true; 
  int* others = new int[k3DNeighbors];
	// Get Neighbors
  GetNeighbourCells(Cell,others);
  // If There Are Zero then False
  for(int loopA=0;loopA<k3DNeighbors;loopA++){
		isInside = ((isInside)&&(others[loopA]>-1));
	} 
  // Deallocate
  delete [] others;
  // RETURN
	return isInside;
}

// GET NEIGHBORS OF A CELL
void MRIScan::GetNeighbourCells(int CurrentCell, int* &cellNeighbors){
  int* coords = new int[3];
  //Get The Coordinates of the Current Cell
  MapIndexToCoords(CurrentCell,coords);
  // Get Neighbor
	// coords[0]
  if ((coords[0]-1)>=0) cellNeighbors[0] = MapCoordsToIndex(coords[0]-1,coords[1],coords[2]);
	else cellNeighbors[0] = -1;	
	// coords[1]
  if((coords[0]+1)<cellTotals[0]) cellNeighbors[1] = MapCoordsToIndex(coords[0]+1,coords[1],coords[2]);
  else cellNeighbors[1] = -1;
  // coords[2]
  if((coords[1]-1)>=0) cellNeighbors[2] = MapCoordsToIndex(coords[0],coords[1]-1,coords[2]);
	else cellNeighbors[2] = -1;
	// coords[3]
  if((coords[1]+1)<cellTotals[1]) cellNeighbors[3] = MapCoordsToIndex(coords[0],coords[1]+1,coords[2]);
	else cellNeighbors[3] = -1;
	// coords[4]
  if((coords[2]-1)>=0) cellNeighbors[4] = MapCoordsToIndex(coords[0],coords[1],coords[2]-1);
  else cellNeighbors[4] = -1;
	// coords[5]
  if((coords[2]+1)<cellTotals[2]) cellNeighbors[5] = MapCoordsToIndex(coords[0],coords[1],coords[2]+1);
  else cellNeighbors[5] = -1;
  // Deallocate
  delete [] coords;
}

// Get Face Starting From Cell and Unit Vector
int MRIScan::GetFacewithCellVector(int CurrentCell, double* UnitVector){
  double tolerance = 5.0e-3;
	int AdjType = 0;
	int resultFace = -1;
	// Check Direction Vector
  if ((fabs(UnitVector[0]-1.0)<tolerance)&&(fabs(UnitVector[1])<tolerance)&&(fabs(UnitVector[2])<tolerance)){
		AdjType = kfacePlusX;
	}else if((fabs(UnitVector[0]+1.0)<tolerance)&&(fabs(UnitVector[1])<tolerance)&&(fabs(UnitVector[2])<tolerance)){
		AdjType = kfaceMinusX;
	}else if((fabs(UnitVector[0])<tolerance)&&(fabs(UnitVector[1]-1.0)<tolerance)&&(fabs(UnitVector[2])<tolerance)){
		AdjType = kfacePlusY;
	}else if((fabs(UnitVector[0])<tolerance)&&(fabs(UnitVector[1]+1.0)<tolerance)&&(fabs(UnitVector[2])<tolerance)){
		AdjType = kfaceMinusY;
	}else if((fabs(UnitVector[0])<tolerance)&&(fabs(UnitVector[1])<tolerance)&&(fabs(UnitVector[2]-1.0)<tolerance)){
		AdjType = kfacePlusZ;
	}else if((fabs(UnitVector[0])<tolerance)&&(fabs(UnitVector[1])<tolerance)&&(fabs(UnitVector[2]+1.0)<tolerance)){
		AdjType = kfaceMinusZ;
	}else{
		throw new MRIMeshCompatibilityException("Internal Error: Problems in GetFacewithCellVector");
	} 
  // Get Adj Face
  resultFace = GetAdjacentFace(CurrentCell,AdjType);
  return resultFace;
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

// Get Global Coords From Local Ones
void MRIScan::GetGlobalCoords(int DimNumber, int SliceNumber, double FaceCoord1, double FaceCoord2, double* &globalCoords){
  switch(DimNumber){
		case 0:
		  // X Direction
      globalCoords[0] = domainSizeMin[0] + cellLength[0]*(SliceNumber);
      globalCoords[1] = domainSizeMin[1] + FaceCoord1;
      globalCoords[2] = domainSizeMin[2] + FaceCoord2;
		  break;
    case 1:
      // Y Direction
      globalCoords[0] = domainSizeMin[0] + FaceCoord1;
      globalCoords[1] = domainSizeMin[1] + cellLength[1]*(SliceNumber);
      globalCoords[2] = domainSizeMin[2] + FaceCoord2;
			break;
    case 2:
      // Z Direction
      globalCoords[0] = domainSizeMin[0] + FaceCoord1;
      globalCoords[1] = domainSizeMin[1] + FaceCoord2;
      globalCoords[2] = domainSizeMin[2] + cellLength[2]*(SliceNumber);
			break;
	}
}

// Transform From Local to Global Face Number
int MRIScan::FaceLocaltoGlobal(int LocalFace, int DimNumber, int SliceNumber){
	int cells1 = 0;
    double cellLength1 = 0.0;
    double cellLength2 = 0.0;
  // If Face is Zero then Return Zero
  if(LocalFace == -1) return -1;
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
}

// Map To Cells Coords
void MRIScan::MapIndexToCoords(int index, int* intCoords){
  int CurrentIndex = index;
  intCoords[2] = (int)(CurrentIndex/(cellTotals[0]*cellTotals[1]));
  CurrentIndex = (CurrentIndex-intCoords[2]*cellTotals[0]*cellTotals[1]);
  intCoords[1] = (int)(CurrentIndex/cellTotals[0]);
  CurrentIndex = CurrentIndex-intCoords[1]*cellTotals[0];
  intCoords[0] = CurrentIndex;
}

// Map From Cells Coords
int MRIScan::MapCoordsToIndex(int i, int j, int k){
	// C++ INDEXES ZERO BASED
  return k*(cellTotals[0]*cellTotals[1])+j*(cellTotals[0])+i;
}

// Map Cell Number
int MRIScan::GetCellNumber(MRIReal* coords){
  // Check Indexes
  int i = MRIUtils::round((coords[0]-domainSizeMin[0])/cellLength[0]);
  int j = MRIUtils::round((coords[1]-domainSizeMin[1])/cellLength[1]);
  int k = MRIUtils::round((coords[2]-domainSizeMin[2])/cellLength[2]);
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
void MRIScan::GetGlobalPermutation(int* &GlobalPerm){
  MRIReal Norm = 0.0;
  MRIReal VelNorm = 0.0;
  int invalidCount = 0;
  // Fill Permutation Vector
  int NumberOfZeroNorms = 0;
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    Norm = MRIUtils::Do3DEucNorm(cellPoints[loopA].position);
    VelNorm = MRIUtils::Do3DEucNorm(cellPoints[loopA].velocity);
    if (Norm<kMathZero) invalidCount++;
    if ((Norm<kMathZero)&&(VelNorm<kMathZero)&&(invalidCount>1)){
      GlobalPerm[loopA] = -1;
      NumberOfZeroNorms++;
    }else{
      GlobalPerm[loopA] = GetCellNumber(cellPoints[loopA].position);
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
int MRIScan::GetAdjacentFace(int globalNodeNumber /*Already Ordered Globally x-y-z*/, int AdjType){
  // Get The Z Coord
  double currentZCoord = cellPoints[globalNodeNumber].position[2]-domainSizeMin[2];
  
	// Find The Node Number in The Current Plane
  int ZCompleteLevels = MRIUtils::round(currentZCoord/cellLength[2]);
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

// Get Local Adjacent Plane
void GetLocalStarFaces(int starNum, int cellsX, int cellsY, 
                       int &bottomFace, int &topFace, int &leftFace, int &rightFace){
  // Find Local Face Number
  bottomFace = (((int)(starNum)/(int)(cellsX+1))-1)*(2*cellsX+1) + ((starNum) % (cellsX+1)) + cellsX + 1;
  topFace    = (((int)(starNum)/(int)(cellsX+1)))  *(2*cellsX+1) + ((starNum) % (cellsX+1)) + cellsX + 1;
  leftFace   = (((int)(starNum)/(int)(cellsX+1)))  *(2*cellsX+1) + ((starNum) % (cellsX+1));
  rightFace  = leftFace+1;

  // Set To Zero the Null Faces
  int iCoord = ((int)(starNum)/(int)(cellsX+1));
  int jCoord = (     (starNum)%     (cellsX+1));
  if(iCoord == 0) bottomFace = -1;
  if(iCoord == cellsY) topFace = -1;
  if(jCoord == 0) leftFace = -1;
  if(jCoord == cellsX) rightFace = -1;
}
