#include <stdlib.h>
#include <math.h>

#include "mriScan.h"
#include "mriStructuredScan.h"
#include "mriConstants.h"
#include "mriUtils.h"
#include "mriThresholdCriteria.h"
#include "mriException.h"
#include "mriExpansion.h"
#include "schMessages.h"

// PRINT FASE ID INDEX
void PrintFaceIDIndexes(std::string fileName, int totalStarFaces, std::vector<int> facesID, std::vector<double> facesCoeffs){
	// Open Output File
	FILE* outFile;
	outFile = fopen(fileName.c_str(),"w");
	// Write Header
  for(int loopA=0;loopA<totalStarFaces;loopA++){
    fprintf(outFile,"%d %d %e\n",loopA,facesID[loopA],facesCoeffs[loopA]);
  }
	// Close Output file
	fclose(outFile);  
}

// Print Residual Vector
void PrintResidualVector(std::string fileName, int totalFaces, double* resVec){
	// Open Output File
	FILE* outFile;
	outFile = fopen(fileName.c_str(),"w");
	// Write Header
  for(int loopA=0;loopA<totalFaces;loopA++){
    fprintf(outFile,"%e\n",resVec[loopA]);
  }
	// Close Output file
	fclose(outFile);
}

// ASSEMBLE STAR SHAPE
void MRIStructuredScan::AssembleConstantPattern(int currentDim, int &totalConstantFaces, std::vector<int> &facesID, std::vector<double> &facesCoeffs){
  // Clear Vectors
  totalConstantFaces = 0;
  facesID.clear();
  facesCoeffs.clear();  
  // Eval The total Number Of Cells
  switch(currentDim){
    case 0: 
      // X
      totalConstantFaces = cellTotals[1]*cellTotals[2]*(cellTotals[0]+1);
      break;
    case 1: 
      // Y
      totalConstantFaces = cellTotals[0]*cellTotals[2]*(cellTotals[1]+1);
      break;
    case 2: 
      // Z
      totalConstantFaces = cellTotals[0]*cellTotals[1]*(cellTotals[2]+1);
      break;
  }
  // Loop On Cells
  int currentComponent = 0;
  int testCount = 0;
	int faceLocation = 0;
	int currentFace = 0;
	bool isLastCell = false;
	int* intCoords = new int[3];
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Get Coordinate of the Present Cell
    MapIndexToCoords(loopA,intCoords);
    // Get Face Type
    switch(currentDim){
      case 0: 
			  faceLocation = kfaceMinusX;
				break;
      case 1: 
			  faceLocation = kfaceMinusY;
				break;
      case 2: 
			  faceLocation = kfaceMinusZ;
				break;
		}
    // Get Current Face
    currentFace = GetAdjacentFace(loopA,faceLocation);
    // Assemble
    currentComponent++;
    facesID.push_back(currentFace);
    facesCoeffs.push_back((1.0/sqrt(totalConstantFaces)));
    // Check If Is Last Cell
    isLastCell = (intCoords[currentDim] == (cellTotals[currentDim]-1));
    if(isLastCell){
      // Get Face Type
      switch(currentDim){
        case 0: 
				  faceLocation = kfacePlusX;
					break;
        case 1: 
				  faceLocation = kfacePlusY;
					break;
        case 2: 
				  faceLocation = kfacePlusZ;
					break;
      }
      // Get Current Face
      currentFace = GetAdjacentFace(loopA,faceLocation);
      // Assemble
      currentComponent++;
      testCount++;
      facesID.push_back(currentFace);
      facesCoeffs.push_back((1.0/sqrt(totalConstantFaces)));
    }
	}
  // Deallocate
  delete [] intCoords;
}

// ASSEMBLE STAR MATRIX
void MRIStructuredScan::AssembleStarMatrix(int &totalFaces, int &totalBasis, double** &starMatrix){
  // Init Rows and Columns
  totalFaces = GetTotalFaces();
  totalBasis = GetTotalBasisNumber();
  // Init Count
  int totalStarFaces = 0;
  std::vector<int> facesID;
  std::vector<double> facesCoeffs;
  int currentRow = 0;
  int currentCount = 0;
  // Get total Number of Vortexes
  int totalVortices = EvalTotalVortex();
  // Allocate and Initialize
  starMatrix = new double*[totalFaces];
  for(int loopA=0;loopA<totalFaces;loopA++) starMatrix[loopA] = new double[totalBasis];
  // Loop On The Three Directions
	int totalSlices = 0;
	int totalStars = 0;
  for(int loopB=0;loopB<totalVortices;loopB++){
    // Find Star Shape
    AssembleStarShape(loopB,totalStarFaces,facesID,facesCoeffs);
    // Fill Matrix
    for(int loopE=0;loopE<totalStarFaces;loopE++){
      currentRow = facesID[loopE];
      starMatrix[currentRow][currentCount] = facesCoeffs[loopE];
    }
    // Update Count
    currentCount++;
  }
}

// GET TOTAL NUMBER OF BASIS
int MRIStructuredScan::GetTotalBasisNumber(){
  int totBasis = 0;
  int totalSlices = 0;
  int totalStars = 0;
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    switch(loopA){
      case 0:
        // YZ Planes
        totalSlices = cellTotals[0];
        totalStars = (cellTotals[1]+1)*(cellTotals[2]+1);
        break;
      case 1:
        // XZ Planes
        totalSlices = cellTotals[1];
        totalStars = (cellTotals[0]+1)*(cellTotals[2]+1);
        break;
      case 2:
        // XY Planes
        totalSlices = cellTotals[2];
        totalStars = (cellTotals[0]+1)*(cellTotals[1]+1);
        break;
    }
    totBasis = totBasis + totalSlices * totalStars;
  }
	return totBasis;
}

// CHECK PERMUTATION
bool CheckPermutation(int size, int* perm){
  bool checkPerm[size];
  for(int loopA=0;loopA<size;loopA++) checkPerm[size] = false;
  // Sort Permutation
  for(int loopA=0;loopA<size;loopA++) checkPerm[perm[loopA]] = true;
  for(int loopA=0;loopA<size;loopA++){
    if (!checkPerm[loopA]) return false;
  }
	return true;
}

// EXPAND STAR SHAPE TO FULL VECTOR
void MRIStructuredScan::ExpandStarShape(int totalStarFaces, int* facesID, double* facesCoeffs, double* &fullStarVector){
  // Get Total Number Of Faces
  int totalFaces = cellTotals[0]*cellTotals[1]*(cellTotals[2]+1)+
                   cellTotals[1]*cellTotals[2]*(cellTotals[0]+1)+
                   cellTotals[2]*cellTotals[0]*(cellTotals[1]+1);
  // Allocate
  fullStarVector = new double[totalFaces];
  for(int loopA=0;loopA<totalFaces;loopA++) fullStarVector[loopA] = 0.0;
  // Fill Array
  int currentFace = 0;
  for(int loopA=0;loopA<totalFaces;loopA++){
    currentFace = facesID[loopA];
    fullStarVector[currentFace] = facesCoeffs[loopA];
  }
}

// EVAL DIVERGENCE
double MRIStructuredScan::EvalMaxDivergence(double* filteredVec){
  double maxDivergence = 0.0;
	double currentDiv = 0.0;
	int faceLocMinus,faceLocPlus;
	int currentFacePlus,currentFaceMinus;
	double fluxPlus,fluxMinus;
  for(int loopA=0;loopA<totalCellPoints;loopA++){
		currentDiv = 0.0;
    for(int loopB=0;loopB<3;loopB++){
      switch(loopB){
        case 0: 
          faceLocPlus = kfacePlusX;
          faceLocMinus = kfaceMinusX;
          break;
        case 1:
          faceLocPlus = kfacePlusY;
          faceLocMinus = kfaceMinusY;
          break;
        case 2:
          faceLocPlus = kfacePlusZ;
          faceLocMinus = kfaceMinusZ;
          break;
      }
      // Get Current Face
      currentFacePlus = GetAdjacentFace(loopA,faceLocPlus);
      currentFaceMinus = GetAdjacentFace(loopA,faceLocMinus);
      // Compute the Two Signs
      fluxPlus = (filteredVec[currentFacePlus]);
      fluxMinus = (filteredVec[currentFaceMinus]);
      // Update
      currentDiv = currentDiv - fluxPlus + fluxMinus;
    }
    // Store Max Absolute Value
    if (fabs(currentDiv)>maxDivergence) maxDivergence = fabs(currentDiv);
  }
  return maxDivergence;
}

// REORDER CELLS
void MRIStructuredScan::ReorderCells(int* Perm){
  // Allocate a copy of the Scan: Check if necessary!!!
  std::vector<MRICell> tempCellPoints;
  tempCellPoints.reserve(totalCellPoints);
  // Initialize
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Conc
    tempCellPoints[loopA].concentration = 0.0;
    // Vel
    tempCellPoints[loopA].velocity[0] = 0.0;
    tempCellPoints[loopA].velocity[1] = 0.0;
    tempCellPoints[loopA].velocity[2] = 0.0;
    // Pos
    tempCellPoints[loopA].position[0] = 0.0;
    tempCellPoints[loopA].position[1] = 0.0;
    tempCellPoints[loopA].position[2] = 0.0;
  }
  // Transfer
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    if(Perm[loopA]>-1){
      // Conc
      tempCellPoints[Perm[loopA]].concentration = cellPoints[loopA].concentration;
      // Vel
      tempCellPoints[Perm[loopA]].velocity[0] = cellPoints[loopA].velocity[0];
      tempCellPoints[Perm[loopA]].velocity[1] = cellPoints[loopA].velocity[1];
      tempCellPoints[Perm[loopA]].velocity[2] = cellPoints[loopA].velocity[2];
      // Pos
      tempCellPoints[Perm[loopA]].position[0] = cellPoints[loopA].position[0];
      tempCellPoints[Perm[loopA]].position[1] = cellPoints[loopA].position[1];
      tempCellPoints[Perm[loopA]].position[2] = cellPoints[loopA].position[2];
    }
  }
  // Copy Back
  double PosNorm;
  double VelNorm;
  int invalidCount = 0;
  int* intCoords = new int[3];
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Conc
    cellPoints[loopA].concentration = tempCellPoints[loopA].concentration;
    // Vel
    cellPoints[loopA].velocity[0] = tempCellPoints[loopA].velocity[0];
    cellPoints[loopA].velocity[1] = tempCellPoints[loopA].velocity[1];
    cellPoints[loopA].velocity[2] = tempCellPoints[loopA].velocity[2];
    // Pos
    PosNorm = MRIUtils::Do3DEucNorm(tempCellPoints[loopA].position);
    VelNorm = MRIUtils::Do3DEucNorm(tempCellPoints[loopA].velocity);
    if(PosNorm<kMathZero) invalidCount++;
    if((PosNorm<kMathZero)&&(VelNorm<kMathZero)&&(invalidCount>1)){
      // Get Coords
      MapIndexToCoords(loopA,intCoords);
      // Assign
      cellPoints[loopA].position[0] = intCoords[0]*cellLength[0] + domainSizeMin[0];
      cellPoints[loopA].position[1] = intCoords[1]*cellLength[1] + domainSizeMin[1];
      cellPoints[loopA].position[2] = intCoords[2]*cellLength[2] + domainSizeMin[2];
    }else{
      cellPoints[loopA].position[0] = tempCellPoints[loopA].position[0];
      cellPoints[loopA].position[1] = tempCellPoints[loopA].position[1];
      cellPoints[loopA].position[2] = tempCellPoints[loopA].position[2];
    }
  }
  // Deallocate
  delete [] intCoords;
}

// RECOVER VELOCITIES FROM FACE FLUXES
void MRIStructuredScan::RecoverCellVelocitiesRT0(bool useBCFilter, double* filteredVec){
  // Variables
  int faceLocPlus = 0;
  double faceArea = 0.0;
  int faceLocMinus = 0;
  int currentFacePlus = 0;
  int currentFaceMinus = 0;
  double avVelocityPlus = 0.0;
  double avVelocityMinus = 0.0;
  // Loop To Assemble Residual Vector
  for(int loopA=0;loopA<totalCellPoints;loopA++){
		// Loop On Faces
		for(int loopB=0;loopB<3;loopB++){
			switch(loopB){
				case 0:
					faceLocPlus = kfacePlusX;
					faceArea = cellLength[1]*cellLength[2];
					faceLocMinus = kfaceMinusX;
					break;
				case 1:
					faceLocPlus = kfacePlusY;
					faceArea = cellLength[0]*cellLength[2];
					faceLocMinus = kfaceMinusY;
					break;
				case 2:
					faceLocPlus = kfacePlusZ;
					faceArea = cellLength[0]*cellLength[1];
					faceLocMinus = kfaceMinusZ;
					break;
			}
			// Get Current Face
			currentFacePlus = GetAdjacentFace(loopA,faceLocPlus);
			currentFaceMinus = GetAdjacentFace(loopA,faceLocMinus);
			// Eval Average Velocity
			avVelocityPlus = (filteredVec[currentFacePlus]/faceArea);
			avVelocityMinus = (filteredVec[currentFaceMinus]/faceArea);
			// Set Correction
			if(useBCFilter){
                cellPoints[loopA].auxVector[loopB] = cellPoints[loopA].velocity[loopB]-0.5*(avVelocityPlus + avVelocityMinus);
			}else{
                cellPoints[loopA].auxVector[loopB] = 0.5*(avVelocityPlus + avVelocityMinus);
			}
	  }
	}
}

// =======================================================
// GET DIMENSION, SLICE AND STAR NUMBER FROM VORTEX NUMBER
// =======================================================
void MRIStructuredScan::getDimensionSliceStarFromVortex(int vortexNumber,int &dimNumber,int &sliceNumber,int &starNumber){
  // Declare
  int totalSlices[3] = {0.0};
  int totalStars[3] = {0.0};
  int vortexOffset = 0;
  
  // Get Dimension,Slice and Star from the cardinality
  // X Direction
  totalSlices[0] = cellTotals[0];
  totalStars[0] = (cellTotals[1]+1)*(cellTotals[2]+1);
  // Y Direction
  totalSlices[1] = cellTotals[1];
  totalStars[1] = (cellTotals[0]+1)*(cellTotals[2]+1);
  // Z Direction
  totalSlices[2] = cellTotals[2];
  totalStars[2] = (cellTotals[0]+1)*(cellTotals[1]+1);

  // Get the total Vortexes for every dimension
  double totalVortexDim1 = totalSlices[0]*totalStars[0];
  double totalVortexDim2 = totalVortexDim1 + totalSlices[1]*totalStars[1];
  double totalVortexDim3 = totalVortexDim2 + totalSlices[2]*totalStars[2];

  // Check the dimension
  if(vortexNumber<totalVortexDim1){
    dimNumber = 0;
    vortexOffset = 0;
  }else if(vortexNumber<totalVortexDim2){
    dimNumber = 1;
    vortexOffset = totalVortexDim1;
  }else if(vortexNumber<totalVortexDim3){
    dimNumber = 2;
    vortexOffset = totalVortexDim2;
  }else{
    throw MRIException("Total Number of Vortices Exceeded.\n");
  }

  // Get Slice and Star Number
  // CHECK !!!
  sliceNumber = (vortexNumber - vortexOffset)/(int)totalSlices[dimNumber];
  starNumber = (vortexNumber - vortexOffset) % (int)totalSlices[dimNumber];

}

// ====================
// ASSEMBLE STAR SHAPES
// ====================
void MRIStructuredScan::AssembleStarShape(int vortexNumber, int &totalFaces, std::vector<int> &facesID, std::vector<double> &facesCoeffs){

  int cells1 = 0;
  int cells2 = 0;
  int localBottomFace = 0;
  int localTopFace = 0;
  int localLeftFace = 0;
  int localRightFace = 0;
  // Faces
  int bottomFace = 0;
  int topFace = 0;
  int leftFace = 0;
  int rightFace = 0;
  // Clear Array and Reset Counters
  totalFaces = 0;
  facesID.clear();
  facesCoeffs.clear();
  int dimNumber = 0;
  int sliceNumber = 0;
  int starNumber = 0;

  // Get Dimension, Slice and Star
  getDimensionSliceStarFromVortex(vortexNumber,dimNumber,sliceNumber,starNumber);

  // Assemble
  switch(dimNumber){
    case 0:
      // YZ
      cells1 = cellTotals[1];
      cells2 = cellTotals[2];
      break;
    case 1:
      // XZ
      cells1 = cellTotals[0];
      cells2 = cellTotals[2];
      break;
    case 2:
      // XY
      cells1 = cellTotals[0];
      cells2 = cellTotals[1];
      break;
  }
  // Get Local Adjacent Plane
  GetLocalStarFaces(starNumber,cells1,cells2,localBottomFace,localTopFace,localLeftFace,localRightFace);
  // Convert To Global Faces
  bottomFace = FaceLocaltoGlobal(localBottomFace,dimNumber,sliceNumber);
  topFace =    FaceLocaltoGlobal(localTopFace,dimNumber,sliceNumber);
  leftFace =   FaceLocaltoGlobal(localLeftFace,dimNumber,sliceNumber);
  rightFace =  FaceLocaltoGlobal(localRightFace,dimNumber,sliceNumber);
  // Insert In Lists
  if(bottomFace>-1){
    MRIUtils::InsertInIntList(bottomFace,totalFaces,facesID);
    facesCoeffs.push_back(-1.0);
  }
  if(topFace>-1){
    MRIUtils::InsertInIntList(topFace,totalFaces,facesID);
    facesCoeffs.push_back(1.0);
  }
  if(leftFace>-1){
    MRIUtils::InsertInIntList(leftFace,totalFaces,facesID);
    facesCoeffs.push_back(1.0);
  }
  if(rightFace>-1){
    MRIUtils::InsertInIntList(rightFace,totalFaces,facesID);
    facesCoeffs.push_back(-1.0);
  }
  // Normalize
  double norm = 0.0;
  for(int loopA=0;loopA<totalFaces;loopA++){
    norm = norm + (facesCoeffs[loopA]*facesCoeffs[loopA]);
  }
  norm = sqrt(norm);
  for(int loopA=0;loopA<totalFaces;loopA++){
    facesCoeffs[loopA] = (facesCoeffs[loopA]/norm);
  }
}

// UPDATE THE RESIDUAL AND RESIDUAL NORM
void UpdateResidualAndFilter(double corrCoeff, int totalFaces, std::vector<int> facesID, std::vector<double> facesCoeffs,
                             double* resVec, double* filteredVels, double &resNorm){
  resNorm = (resNorm*resNorm);
  int currentFaceID = 0;
  double currentFaceCoeff = 0.0;
  for(int loopA=0;loopA<totalFaces;loopA++){
    currentFaceID = facesID[loopA];
    currentFaceCoeff = facesCoeffs[loopA];
    // Update Residual Norm
    resNorm = resNorm - (resVec[currentFaceID]*resVec[currentFaceID]) + 
                        (resVec[currentFaceID]-corrCoeff*currentFaceCoeff)*(resVec[currentFaceID]-corrCoeff*currentFaceCoeff);
    // UPDATE RESIDUAL AND RECONSTRUCTED FIELD
    // Update Residual
    resVec[currentFaceID] = resVec[currentFaceID] - corrCoeff * currentFaceCoeff;
    // Update Filtered Vector
    filteredVels[currentFaceID] = filteredVels[currentFaceID] + corrCoeff * currentFaceCoeff;
  }
  resNorm = sqrt(fabs(resNorm));
}

// EVAL CORRELATION COEFFICIENT
double EvalCoerrelationCoefficient(double* resVec, int totalFaces, const std::vector<int> &facesID, const std::vector<double> &facesCoeffs){
  double corrCoeff = 0.0;
  int currentFaceID = 0;
  double currentCoeff = 0.0;
  for(int loopA=0;loopA<totalFaces;loopA++){
    currentFaceID = facesID[loopA];
    currentCoeff = facesCoeffs[loopA];
    corrCoeff += resVec[currentFaceID] * currentCoeff;
  }
  return corrCoeff;
}

// GET A STRING ASSOCIATED TO THE FACE LOCATION
std::string GetFaceLocationString(int faceLocation){
  switch(faceLocation){
    case kfacePlusX:  
      return std::string("PlusX");
      break;
    case kfaceMinusX: 
      return std::string("MinusX");
      break;
    case kfacePlusY:  
      return std::string("PlusY");
      break;
    case kfaceMinusY: 
      return std::string("MinusY");
      break;
    case kfacePlusZ:  
      return std::string("PlusZ");
      break;
    case kfaceMinusZ: 
      return std::string("MinusZ");
      break;
  }
  return std::string("");
}

// Assemble Face Flux Vectors
void MRIStructuredScan::AssembleResidualVector(bool useBCFilter, MRIThresholdCriteria thresholdCriteria,
                                                 int &totalFaces, double* &resVec, double* &filteredVec, double &resNorm){
  bool   continueToProcess = false;
  double currentValue = 0.0;
  double currentVelX = 0.0;
  double currentVelY = 0.0;
  double currentVelZ = 0.0;
  int    faceLocation = 0;
  double faceArea = 0.0;
  double faceComponent = 0.0;
  int    currentFace = 0;
  bool   checkPassed = false;
  // Get Total Number Of Faces
  totalFaces = cellTotals[0]*cellTotals[1]*(cellTotals[2]+1)+
               cellTotals[1]*cellTotals[2]*(cellTotals[0]+1)+
               cellTotals[2]*cellTotals[0]*(cellTotals[1]+1);
  // Allocate
  resVec = new double[totalFaces];
  filteredVec = new double[totalFaces];
  int* resID = new int[totalFaces];
  // Initialize
  for(int loopA=0;loopA<totalFaces;loopA++){
    resVec[loopA] = 0.0;
  }
  for(int loopA=0;loopA<totalFaces;loopA++){
    resID[loopA] = 0;
  }
  for(int loopA=0;loopA<totalFaces;loopA++){
    filteredVec[loopA] = 0.0;
  }
  // Loop To Assemble Residual Vector
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Check for BC
    if(useBCFilter){
      currentValue = cellPoints[loopA].getQuantity(thresholdCriteria.thresholdQty);
      continueToProcess = thresholdCriteria.MeetsCriteria(currentValue);
    }else{
      continueToProcess = true;
    }
    if(continueToProcess){
      // Store Local Values
      currentVelX = cellPoints[loopA].velocity[0];
      currentVelY = cellPoints[loopA].velocity[1];
      currentVelZ = cellPoints[loopA].velocity[2];
      // Loop On Faces
      for(int loopB=0;loopB<k3DNeighbors;loopB++){
        switch(loopB){
          case 0:
            faceLocation = kfacePlusX;
            faceArea = cellLength[1]*cellLength[2];
            faceComponent = currentVelX;
            break;
          case 1:
            faceLocation = kfaceMinusX;
            faceArea = cellLength[1]*cellLength[2];
            faceComponent = currentVelX;
            break;
          case 2:
            faceLocation = kfacePlusY;
            faceArea = cellLength[0]*cellLength[2];
            faceComponent = currentVelY;
            break;
          case 3:
            faceLocation = kfaceMinusY;
            faceArea = cellLength[0]*cellLength[2];
            faceComponent = currentVelY;
            break;
          case 4:
            faceLocation = kfacePlusZ;
            faceArea = cellLength[0]*cellLength[1];
            faceComponent = currentVelZ;
            break;
          case 5:
            faceLocation = kfaceMinusZ;
            faceArea = cellLength[0]*cellLength[1];
            faceComponent = currentVelZ;
            break;
        }        
        // Get Current Face
        currentFace = GetAdjacentFace(loopA,faceLocation);
        // Assemble
        resVec[currentFace] = resVec[currentFace] + faceArea * faceComponent;
        resID[currentFace]++;
      }
    }
  }
  // Check Faces
  for(int loopA=0;loopA<totalFaces;loopA++){
    if(useBCFilter) checkPassed = (resID[loopA]>2);
    else checkPassed = (resID[loopA]<1)||(resID[loopA]>2);
    if(checkPassed){
      std::string currentMsgs = "Internal: Wrong Face Connectivity, Face: " + MRIUtils::IntToStr(loopA)+ "; Connectivity: " + MRIUtils::IntToStr(resID[loopA])+".";
      throw new MRIMeshCompatibilityException(currentMsgs.c_str());
    }
  }
  /*// Count Faces
  int SingleConnFaces = 0;
  int DoubleConnFaces = 0;
  for(int loopA=0;loopA<totalFaces;loopA++){
    if (resID[loopA] == 1){
      SingleConnFaces++;
    }else{
      DoubleConnFaces++;
    }
  }  
  printf("Single Count: %d; Double Count %d\n",SingleConnFaces,DoubleConnFaces);*/
  
  // Divide By the Number Of Faces
  for(int loopA=0;loopA<totalFaces;loopA++){
    //printf("Face %d; Connected %d\n",loopA,resID[loopA]);
    //printf("Face %d; Flux %e\n",loopA,resVec[loopA]);
    if(resID[loopA]>0) resVec[loopA] = ((double)resVec[loopA]/(double)resID[loopA]);
    else resVec[loopA] = 0.0;
  }
  // Find Initial Residual Norm
  resNorm = 0.0;
  for(int loopA=0;loopA<totalFaces;loopA++){
    resNorm = resNorm + (resVec[loopA]*resVec[loopA]);
  }
  resNorm = sqrt(resNorm);
  // Deallocate
  delete [] resID;
}

// =================
// PHYSICS FILTERING
// =================
void MRIScan::applySMPFilter(MRIOptions* options){
  // Initialization
  int totalFaces = 0;
  double* resVec = NULL;
  double* filteredVec = NULL;
  std::vector<int> facesID;
  std::vector<double> facesCoeffs;
  double corrCoeff = 0.0;
  int totalSlices = 0;
  int totalStars = 0;
  int totalStarFaces = 0;
  double maxDivergence = 0.0;
  double maxNormError = 0.0;
  double maxAngleError = 0.0;

  // Set up Norms
  double resNorm = 0.0;
  double relResNorm = 0.0;
  double twoNorm = 0.0;
  double relTwoNorm = 0.0;
  
  // Assemble Face Flux Vectors
  AssembleResidualVector(options->useBCFilter,options->thresholdCriteria,totalFaces,resVec,filteredVec,resNorm);
  
  // Print Residual Vector
  //PrintResidualVector("resVectorFile.txt",totalFaces,resVec);

  // Initial Residual
  WriteSchMessage("\n");
  if (options->useBCFilter){
    WriteSchMessage("FILTER ALGORITHM - BC - Step: "+MRIUtils::FloatToStr(scanTime)+" ---------------------------\n");
  }else{
    WriteSchMessage("FILTER ALGORITHM - FULL - Step "+MRIUtils::FloatToStr(scanTime)+" ---------------------------\n");
  }

  // START CLOCK
  const clock_t begin_time = clock();
  
  WriteSchMessage("Initial Residual Norm: "+MRIUtils::FloatToStr(resNorm)+"\n");

  // Initialize Expansion
  MRIExpansion* bcExpansion = NULL;
  int totalVortexes = EvalTotalVortex();
  if(!options->useBCFilter){
    expansion = new MRIExpansion(totalVortexes);
  }else{
    bcExpansion = new MRIExpansion(totalVortexes);
  }

  // Apply MP Filter
  bool converged = false;
  int itCount = 0;
  double oldResNorm = resNorm;
  double oldTwoNorm = twoNorm;
  // TEMPORARY!!!
  bool stopIterations = false;
  // Set Tolerance
  double itTol = options->tolerance;
  // Initialize Component Count
  int componentCount = 0;
  // Start Filter Loop
  while((!converged)&&(itCount<options->maxIterations)&&(!stopIterations)){
    // Update Iteration Count
    itCount++;
    // Loop On The Three Directions
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      // Correlate Constant Patterns
      if(options->useConstantPatterns){
        // Find Star Shape
        AssembleConstantPattern(loopB/*Current Dimension*/,totalStarFaces,facesID,facesCoeffs);
        // PRINT FASE ID INDEX
        //PrintFaceIDIndexes("ConstantFaceIDs.txt",totalStarFaces,facesID,facesCoeffs);
        // Find Correlation
        corrCoeff = EvalCoerrelationCoefficient(resVec,totalStarFaces,facesID,facesCoeffs);
        // Store Correlation coefficient in Expansion
        // sistemare per BCFilter !!!
        if(!options->useBCFilter){
          expansion->constantFluxCoeff[loopB] += corrCoeff;
        }else{
          bcExpansion->constantFluxCoeff[loopB] += corrCoeff;
        }
        // Update Residual
        UpdateResidualAndFilter(corrCoeff,totalStarFaces,facesID,facesCoeffs,resVec,filteredVec,resNorm);
      }
    }
    // Init componentCount
    componentCount = -1;
    for(int loopB=0;loopB<totalVortexes;loopB++){
      // Increment the current component
      componentCount++;
      // Find Star Shape
      AssembleStarShape(loopB,totalStarFaces,facesID,facesCoeffs);
      // Find Correlation
      corrCoeff = EvalCoerrelationCoefficient(resVec,totalStarFaces,facesID,facesCoeffs);
      // Store Correlation coefficient in Expansion
      if(!options->useBCFilter){
        expansion->vortexCoeff[componentCount] += corrCoeff;
      }else{
        bcExpansion->vortexCoeff[componentCount] += corrCoeff;
      }
      // Update Residual
      UpdateResidualAndFilter(corrCoeff,totalStarFaces,facesID,facesCoeffs,resVec,filteredVec,resNorm);
    }

    // Eval Two-Norm of the Coefficient Vector
    if(!options->useBCFilter){
      twoNorm = expansion->Get2Norm(false);
    }else{
      twoNorm = bcExpansion->Get2Norm(false);
    }

    // Eval Relative Residual Norm
    if(fabs(oldResNorm)>kMathZero){
      relResNorm = fabs((resNorm-oldResNorm)/(oldResNorm));
    }else{
      relResNorm = 0.0;
    }

    // Eval Relative Coefficient Two-Norm
    if(fabs(oldTwoNorm)>kMathZero){
      relTwoNorm = fabs((twoNorm-oldTwoNorm)/(oldTwoNorm));
    }else{
      relTwoNorm = 0.0;
    }

    // WRITE MESSAGE AT EVERY INTERATION
    WriteSchMessage("It: " + MRIUtils::IntToStr(itCount) + "; ABS Res: "+MRIUtils::FloatToStr(resNorm)+"; Rel: " + MRIUtils::FloatToStr(relResNorm) +
                    "; Coeff 2-Norm: "+MRIUtils::FloatToStr(twoNorm)+"; Rel 2-Norm: " + MRIUtils::FloatToStr(relTwoNorm)+"\n");

    // Check Convergence
    //converged = (itCount>1);
    if(itCount>1){
      if(oldResNorm<kMathZero) converged = true;
      else converged = (fabs((resNorm-oldResNorm)/(oldResNorm))<itTol);
    }else converged = false;

    // Update Norm
    oldResNorm = resNorm;
    oldTwoNorm = twoNorm;
  }

  // Iterations Stopped By User
  // TODO: ALLOW THE USER TO STOP THE ITERATIONS!!!
  if(stopIterations){
    WriteSchMessage("ITERATIONS STOPPED BY USER.");
  }

  // Final Residual
  WriteSchMessage("Final Residual Norm: " + MRIUtils::FloatToStr(resNorm) + "\n");

  // Check If the Flux Is Locally Conservative
  maxDivergence = EvalMaxDivergence(filteredVec);

  // Write Divergence Message
  WriteSchMessage("Max Divergence: " + MRIUtils::FloatToStr(maxDivergence) + "\n");

  // Make Diffence between Coefficient Expansions
  if(options->useBCFilter){
    for(int loopA=0;loopA<3;loopA++){
      expansion->constantFluxCoeff[loopA] -= bcExpansion->constantFluxCoeff[loopA];
    }
    for(int loopA=0;loopA<totalVortexes;loopA++){
      expansion->vortexCoeff[loopA] -= bcExpansion->vortexCoeff[loopA];
    }
  }

  // Recover Velocities from Face Fluxes
  RecoverCellVelocitiesRT0(options->useBCFilter,filteredVec);

  // WRITE CPU TIME AND NUMBER OF ITERATIONS
  float totalCPUTime = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
  WriteSchMessage("Total Iterations " + MRIUtils::IntToStr(itCount) + "; Total CPU Time: " + MRIUtils::FloatToStr(totalCPUTime) + "\n");

  // Eval Magniture and Angle Error
  RecoverGlobalErrorEstimates(maxNormError,maxAngleError);

  // Write Final Residual
  WriteSchMessage("Average Magnitude Error: " + MRIUtils::FloatToStr(maxNormError) + "\n");
  WriteSchMessage("Average Angular Error: " + MRIUtils::FloatToStr(maxAngleError) + "\n");

  // Close Filter Comments
  WriteSchMessage("------------------------------------------------------\n");

  // Update Velocities
  //UpdateVelocities(GlobalData);

  // Filtered Velocities Are Available
  //AreVelocityFiltered:=TRUE;
  
  // Deallocate
  delete [] resVec;
  delete [] filteredVec;
}

// ======================================
// REBUILD FROM EXPANSION COEFFICIENT
// ======================================
void MRIStructuredScan::RebuildFromExpansion(MRIExpansion* expansion,bool useConstantFlux){
  // Check Size Compatibility

  // RECONSTRUCTION
  int totalStarFaces = 0;
  int currFaceID = 0;
  double currFaceCoeff = 0.0;
  std::vector<int> facesID;
  std::vector<double> facesCoeffs;
  // Get total number of vortexes
  int totalVortex = EvalTotalVortex();
  double* faceFluxVec = new double[totalVortex+3];

  // INITIALIZE
  for(int loopA=0;loopA<totalVortex+3;loopA++){
    faceFluxVec[loopA] = 0.0;
  }

  // GLOBAL ATOMS
  if(useConstantFlux){
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      // Find Star Shape
      AssembleConstantPattern(loopB/*Current Dimension*/,totalStarFaces,facesID,facesCoeffs);
      // Add to faces
      for(int loopC=0;loopC<totalStarFaces;loopC++){
        currFaceID = facesID[loopC];
        currFaceCoeff = facesCoeffs[loopC];
        faceFluxVec[currFaceID] += currFaceCoeff*expansion->constantFluxCoeff[loopB];
      }
    }
  }

  // VORTEX ATOMS
  int componentCount = -1;
  for(int loopB=0;loopB<expansion->totalVortices;loopB++){
    // Increment the current component
    componentCount++;
    // Find Star Shape
    AssembleStarShape(loopB,totalStarFaces,facesID,facesCoeffs);
    // Find Correlation
    for(int loopE=0;loopE<totalStarFaces;loopE++){
      currFaceID = facesID[loopE];
      currFaceCoeff = facesCoeffs[loopE];
      faceFluxVec[currFaceID] += currFaceCoeff*expansion->vortexCoeff[componentCount];
    }
  }

  // Recover Velocities from Face Fluxes
  RecoverCellVelocitiesRT0(false,faceFluxVec);

  // Update Velocities
  //UpdateVelocities();

  // DELETE ARRAY
  delete [] faceFluxVec;
}

/*{Physics Filtering}
Function PerformPhysicsFilteringLSMR(IsBCFilter: Boolean;
                                     ThresholdCriterion: ThresholdCritRecord;
                                     Var GlobalData: GlobalDataRecord): Integer4;
Const
  NumberOfDimensions = 3;
Var
  LoopA,LoopB,LoopC,LoopD: Integer4;
  ErrorInt: Integer4;
  {Direct Permutation For Ordering}
  DirectPerm: Integer4Array;
  {Inverse Permutation For Ordering}
  InversePerm: Integer4Array;
  {Residual}
  ResVec: DoubleArray;
  ResNorm: Double;
  {Reconstructed Vector}
  FilteredVec: DoubleArray;
  FullStarVector: DoubleArray;
  {Star Patterns}
  TotalStarFaces: Integer4;
  FacesID: Integer4Array;
  FacesCoeffs: DoubleArray;
  {Correlation Coeff}
  CorrCoeff: Double;
  {Total Slices and Stars}
  TotalSlices,TotalStars: Integer4;
  {Error Estimates}
  MaxNormError,MaxAngleError: Double;
  {Max Divergence}
  MaxDivergence: Double;
  IsValidPerm: Boolean;
  CurrentIDX,TotalIDX: Integer4;
  {LSMR Options}
  LSMROptions: LSMROptionRecord;
  ErrorString: String;
  TotalFaces: Integer4;
  TotalBasis: Integer4;
  Converged: Boolean;
  StarMatrix: Double2DArray;
Const
  UseFullLS = FALSE;
Begin
  {Init Result}
  Result:=MRI_IO_NoError;

  {Export Nodes}
  //ExportNodesToFile('G:\Temp\Nodes.csv',GlobalData);

  {Assemble Face Flux Vectors}
  ErrorInt:=AssembleResidualVector(GlobalData,IsBCFilter,ThresholdCriterion,TotalFaces,ResVec,FilteredVec,ResNorm);
  If (ErrorInt<>MRI_IO_NoError) Then
  Begin
    Result:=ErrorInt;
    Exit;
  End;

  {Get the Total Number Of Basis}
  TotalBasis:=GetTotalBasisNumber(GlobalData);

  {Initial Residual}
  WriteMsg('');
  WriteMsg('FILTER ALGORITHM ----------------------------------------');

  {Choose How to Solve}
  If UseFullLS Then
  Begin
    {Use Full LS}
    {Assemble Star Matrix}
    AssembleStarMatrix(GlobalData,TotalFaces,TotalBasis,StarMatrix);
    {Solve With Least Squares}
    ErrorString:=LeastSquaresSolve(TotalFaces,TotalBasis,StarMatrix,ResVec,FilteredVec);
    If (ErrorString<>'') Then
    Begin
      Result:=MRI_LSMR_Problems;
      Exit;
    End;
  End Else Begin
    {Solve With LSMR}
    LSMROptions.MaxIt:=1000;
    LSMROptions.ANormTol:=1.0e-6*ResNorm;
    LSMROptions.BNormTol:=1.0e-6*ResNorm;
    LSMROptions.WriteMsgs:=TRUE;
    LSMROptions.IterationPrintStep:=1;
    ErrorString:=LSMRSolve(LSMROptions,TotalFaces,TotalBasis,ResVec,Converged,FilteredVec);
    If (ErrorString<>'') Then
    Begin
      MessageDlg('Problem: LSMR Solver',mtError,[mbOK],0);
      Result:=MRI_LSMR_Problems;
      Exit;
    End;
  End;

  {Destroy Progress Box}
  DestroyProgressBox;

  {Check If the Flux Is Locally Conservative}
  EvalMaxDivergence(GlobalData,FilteredVec,MaxDivergence);

  {Write Divergence Message}
  WriteMsg('Max Divergence: '+FloatToStr(MaxDivergence));

  {Recover Velocities from Face Fluxes}
  RecoverCellVelocitiesRT0(IsBCFilter,FilteredVec,GlobalData);
  //RecoverCellVelocitiesRT0(ResVec,GlobalData);

  {Eval Magniture and Angle Error}
  RecoverGlobalErrorEstimates(GlobalData,MaxNormError,MaxAngleError);

  {Final Residual}
  WriteMsg('Max Norm Error: '+FloatToStr(MaxNormError));
  WriteMsg('Max Angle Error: '+FloatToStr(MaxAngleError));

  {Close Filter Comments}
  WriteMsg('------------------------------------------------------');
  WriteMsg('');

  {Filtered Velocities Are Available}
  AreVelocityFiltered:=TRUE;
End;


end.*/


