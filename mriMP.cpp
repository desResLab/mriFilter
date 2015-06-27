#include <stdlib.h>
#include <math.h>

#include "mriScan.h"
#include "mriTypes.h"
#include "mriStructuredScan.h"
#include "mriConstants.h"
#include "mriUtils.h"
#include "mriThresholdCriteria.h"
#include "mriException.h"
#include "mriExpansion.h"
#include "schMessages.h"
#include "mriUtils.h"

// ===================
// PRINT FASE ID INDEX
// ===================
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

// =====================
// PRINT RESIDUAL VECTOR
// =====================
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

// ===================
// ASSEMBLE STAR SHAPE
// ===================
void MRIStructuredScan::AssembleConstantPattern(int currentDim, int &totalConstantFaces, std::vector<int> &facesID, std::vector<double> &facesCoeffs){

  // Clear Vectors
  totalConstantFaces = 0;
  facesID.clear();
  facesCoeffs.clear();
  std::vector<int> orientation;

  // Loop Over Faces
  for(size_t loopA=0;loopA<faceConnections.size();loopA++){
    if(fabs(faceNormal[loopA][currentDim])>kMathZero){
      facesID.push_back(loopA);
      if(faceNormal[loopA][currentDim]>kMathZero){
        orientation.push_back(1);
      }else{
        orientation.push_back(-1);
      }

    }
  }

  // Fill Coefficients
  for(size_t loopA=0;loopA<facesID.size();loopA++){
    facesCoeffs.push_back((orientation[loopA]/sqrt((double)facesID.size())));
  }

  // Update Counter
  totalConstantFaces = facesID.size();

}

// ====================
// ASSEMBLE STAR MATRIX
// ====================
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

// =========================
// GET TOTAL NUMBER OF BASIS
// =========================
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

// =================
// CHECK PERMUTATION
// =================
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

// ==============================
// EVAL DIVERGENCE FOR EVERY CELL
// ==============================
double MRIStructuredScan::EvalMaxDivergence(double* filteredVec){
  double maxDivergence = 0.0;
  double currentDiv = 0.0;
  int currFace = 0;
  double extNormal[3] = {0.0};
  double normalSign = 0.0;
  // Loop on cells
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    currentDiv = 0.0;
    // Loop on faces
    for(size_t loopB=0;loopB<cellFaces[loopA].size();loopB++){
      // Get Current Face
      currFace = cellFaces[loopA][loopB];
      // Get External Normal
      getExternalFaceNormal(loopA,loopB,extNormal);
      // Get Sign
      normalSign = 0.0;
      for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
        normalSign += extNormal[loopC] * faceNormal[currFace][loopC];
      }
      round(normalSign);
      // Get Sign
      currentDiv += filteredVec[cellFaces[loopA][loopB]] * normalSign;
    }
    // Store Max Absolute Value
    if (fabs(currentDiv)>maxDivergence){
      maxDivergence = fabs(currentDiv);
    }
  }
  return maxDivergence;
}

// ====================================
// EVAL CELL DIVERGENCE FOR A GIVEN QTY
// ====================================
MRIDoubleVec MRIStructuredScan::evalCellDivergences(MRIDoubleVec faceVec){
  MRIDoubleVec cellDivs;
  double currentDiv = 0.0;
  int currFace = 0;
  double extNormal[3] = {0.0};
  double normalSign = 0.0;

  // Loop on cells
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    currentDiv = 0.0;
    // Loop on faces
    for(size_t loopB=0;loopB<cellFaces[loopA].size();loopB++){
      // Get Current Face
      currFace = cellFaces[loopA][loopB];
      // Get External Normal
      getExternalFaceNormal(loopA,loopB,extNormal);
      // Get Sign
      normalSign = 0.0;
      for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
        normalSign += extNormal[loopC] * faceNormal[currFace][loopC];
      }
      round(normalSign);
      // Get Sign
      currentDiv += faceVec[cellFaces[loopA][loopB]] * normalSign;
    }
    // Store Divergence
    cellDivs.push_back(currentDiv);
  }
  return cellDivs;
}


// =============
// REORDER CELLS
// =============
void MRIStructuredScan::ReorderCells(std::vector<int> Perm){
  // Allocate a copy of the Scan: Check if necessary!!!
  std::vector<MRICell> tempCellPoints;
  tempCellPoints.resize(totalCellPoints);
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
  int intCoords[3] = {0};
  double Pos[3] = {0.0};
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
    if(PosNorm<kMathZero){
      invalidCount++;
    }
    if((PosNorm<kMathZero)&&(VelNorm<kMathZero)&&(invalidCount>1)){
      // Get Coords
      MapIndexToCoords(loopA,intCoords);
      // Get Position
      MapCoordsToPosition(intCoords,true,Pos);
      // Assign
      cellPoints[loopA].position[0] = Pos[0];
      cellPoints[loopA].position[1] = Pos[1];
      cellPoints[loopA].position[2] = Pos[2];
    }else{
      cellPoints[loopA].position[0] = tempCellPoints[loopA].position[0];
      cellPoints[loopA].position[1] = tempCellPoints[loopA].position[1];
      cellPoints[loopA].position[2] = tempCellPoints[loopA].position[2];
    }
  }
}

// ===================================
// RECOVER VELOCITIES FROM FACE FLUXES
// ===================================
void MRIStructuredScan::RecoverCellVelocitiesRT0(bool useBCFilter, double* filteredVec){
  // Variables
  int locFace1 = 0;
  int locFace2 = 0;
  int firstFaceID = 0;
  int secondFaceID = 0;
  double firstAvVelocity = 0.0;
  double secondAvVelocity = 0.0;
  // Loop To Assemble Residual Vector
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      // Get Local Face Numbers
      switch(loopB){
        case 0:
          locFace1 = 2;
          locFace2 = 3;
          break;
        case 1:
          locFace1 = 4;
          locFace2 = 5;
          break;
        case 2:
          locFace1 = 0;
          locFace2 = 1;
          break;
      }
      firstFaceID = cellFaces[loopA][locFace1];
      secondFaceID = cellFaces[loopA][locFace2];
      // Eval Average Velocity
      firstAvVelocity = (filteredVec[firstFaceID] * faceNormal[firstFaceID][loopB]/faceArea[firstFaceID]);
      secondAvVelocity = (filteredVec[secondFaceID]* faceNormal[secondFaceID][loopB]/faceArea[secondFaceID]);
      // Set Correction
      if(useBCFilter){
        cellPoints[loopA].auxVector[loopB] = cellPoints[loopA].velocity[loopB]-0.5*(firstAvVelocity + secondAvVelocity);
      }else{
        cellPoints[loopA].auxVector[loopB] = 0.5*(firstAvVelocity + secondAvVelocity);
      }
    }
  }
}

// =======================================================
// GET DIMENSION, SLICE AND STAR NUMBER FROM VORTEX NUMBER
// =======================================================
void MRIStructuredScan::getDimensionSliceStarFromVortex(int vortexNumber,int &dimNumber,int &sliceNumber,int &starNumber){
  // Declare
  int totalSlices[3] = {0};
  int totalStars[3] = {0};
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
  // Declare
  int currFace = 0;
  double currCoeff = 0.0;

  // Clear Array and Reset Counters
  totalFaces = 0;
  //facesID.reserve(10);
  //facesCoeffs.reserve(10);

  // Loop Over All Faces
  double norm = 0.0;
  for(size_t loopB=0;loopB<edgeFaces[vortexNumber].size();loopB++){
    currFace = abs(edgeFaces[vortexNumber][loopB])-1;
    currCoeff = (double)edgeFaces[vortexNumber][loopB]/(double)fabs(edgeFaces[vortexNumber][loopB]);

    // New Code
    totalFaces++;
    facesID[totalFaces-1] = currFace;
    facesCoeffs[totalFaces-1] = currCoeff;

    norm += currCoeff*currCoeff;

    //MRIUtils::InsertInIntList(currFace,totalFaces,facesID);
    //facesCoeffs.push_back(currCoeff);
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
double EvalCorrelationCoefficient(double* resVec, int totalFaces, const std::vector<int> &facesID, const std::vector<double> &facesCoeffs){
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
void MRIStructuredScan::AssembleResidualVector(bool useBCFilter, MRIThresholdCriteria* thresholdCriteria,
                                                 int &totalFaces, double* &resVec, double* &filteredVec, double &resNorm){
  bool   continueToProcess = false;
  double currentValue = 0.0;
  double faceComponent = 0.0;
  int    currentFace = 0;
  double currFaceArea = 0.0;
  bool   checkPassed = false;

  // Get Total Number Of Faces
  totalFaces = faceConnections.size();

  // Allocate
  resVec = new double[totalFaces];
  filteredVec = new double[totalFaces];
  int* resID = new int[totalFaces];

  // Initialize
  for(int loopA=0;loopA<totalFaces;loopA++){
    resVec[loopA] = 0.0;
    resID[loopA] = 0;
    filteredVec[loopA] = 0.0;
  }

  // Loop To Assemble Residual Vector
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Check for BC
    if(useBCFilter){
      currentValue = cellPoints[loopA].getQuantity(thresholdCriteria->thresholdQty);
      continueToProcess = thresholdCriteria->MeetsCriteria(currentValue);
    }else{
      continueToProcess = true;
    }
    if(continueToProcess){
      // Loop On Faces
      for(int loopB=0;loopB<k3DNeighbors;loopB++){
        // Get Current Face
        currentFace = cellFaces[loopA][loopB];
        // Get Face Area
        currFaceArea = faceArea[currentFace];
        // Get Normal Veclocity
        faceComponent = 0.0;
        for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
          faceComponent += cellPoints[loopA].velocity[loopC] * faceNormal[currentFace][loopC];
        }
        // Assemble
        resVec[currentFace] = resVec[currentFace] + currFaceArea * faceComponent;
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
void MRIStructuredScan::applySMPFilter(MRIOptions* options, bool isBC, MRICommunicator* comm){

  printf("SMPFILTER processor %d\n",comm->currProc);
  // Initialization
  int totalFaces = 0;
  double* resVec = nullptr;
  double* filteredVec = nullptr;
  std::vector<int> facesID;
  std::vector<double> facesCoeffs;
  double corrCoeff = 0.0;
  int totalStarFaces = 0;
  int mpiError = 0;
  double maxDivergence = 0.0;
  double maxNormError = 0.0;
  double maxAngleError = 0.0;

  // Set up Norms
  double resNorm = 0.0;
  double relResNorm = 0.0;
  double twoNorm = 0.0;
  double relTwoNorm = 0.0;

  // Init Time Counters
  float assembleRes_BeginTime = 0.0;
  float assembleRes_TotalTime = 0.0;

  float constPattern_BeginTime = 0.0;
  float constPattern_TotalTime = 0.0;

  float assembleStar_BeginTime = 0.0;
  float assembleStar_TotalTime = 0.0;

  float correlateStar_BeginTime = 0.0;
  float correlateStar_TotalTime = 0.0;

  float updateStar_BeginTime = 0.0;
  float updateStar_TotalTime = 0.0;

  // Determine the Minimum and Maximum Face number of current processor
  int minFaceOnProc = comm->currProc * int((cellFaces.size() - 1)/(comm->totProc + 1));
  int maxFaceOnProc = (comm->currProc + 1) * int((cellFaces.size() - 1)/(comm->totProc + 1))-1;
  if(comm->currProc == (comm->totProc-1)){
    maxFaceOnProc = cellFaces.size();
  }

  // All processes are waiting for the root to read the files
  mpiError = MPI_Barrier(comm->mpiComm);
  MRIUtils::checkMpiError(mpiError);

  printf("Curr Proc: %d, Tot Proc: %d, MinFace: %d, MaxFace: %d\n",comm->currProc,comm->totProc,minFaceOnProc,maxFaceOnProc);

  // Assemble Face Flux Vectors
  assembleRes_BeginTime = clock();
  AssembleResidualVector(isBC,options->thresholdCriteria,totalFaces,resVec,filteredVec,resNorm);
  assembleRes_TotalTime += float( clock () - assembleRes_BeginTime ) /  CLOCKS_PER_SEC;


  // Print Residual Vector
  //PrintResidualVector(std::string("resVectorFile_"+std::to_string(comm->currProc)+".txt").c_str(),totalFaces,resVec);

  // Initial Residual
  WriteSchMessage("\n");
  if (isBC){
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
  if(!isBC){
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
  double itTol = options->itTol;
  // Initialize Component Count
  int componentCount = 0;

  // Start Filter Loop
  while((!converged)&&(itCount<options->maxIt)&&(!stopIterations)){

    // Update Iteration Count
    itCount++;

    constPattern_BeginTime = clock();

    // LOOP ON THE THREE DIRECTIONS
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      // Correlate Constant Patterns
      if(options->useConstantPatterns){
        // Find Star Shape
        AssembleConstantPattern(loopB,totalStarFaces,facesID,facesCoeffs);
        // PRINT FASE ID INDEX
        //PrintFaceIDIndexes("ConstantFaceIDs.txt",totalStarFaces,facesID,facesCoeffs);
        // Find Correlation
        corrCoeff = EvalCorrelationCoefficient(resVec,totalStarFaces,facesID,facesCoeffs);
        // Store Correlation coefficient in Expansion
        // sistemare per BCFilter !!!
        if(!isBC){
          expansion->constantFluxCoeff[loopB] += corrCoeff;
        }else{
          bcExpansion->constantFluxCoeff[loopB] += corrCoeff;
        }
        // Update Residual
        UpdateResidualAndFilter(corrCoeff,totalStarFaces,facesID,facesCoeffs,resVec,filteredVec,resNorm);
      }
    }

    //printf("RESIDUAL NORM AFTER CONSTANT: %e, proc: %d\n",resNorm,comm->currProc);

    constPattern_TotalTime = float( clock () - constPattern_BeginTime ) /  CLOCKS_PER_SEC;

    // Reserve these two
    facesID.resize(10);
    facesCoeffs.resize(10);

    // LOOP ON VORTEXES
    componentCount = -1;
    for(int loopB=0;loopB<totalVortexes;loopB++){
      // Increment the current component
      componentCount++;

      // ASSEMBLE STAR
      assembleStar_BeginTime = clock();
      AssembleStarShape(loopB,totalStarFaces,facesID,facesCoeffs);
      assembleStar_TotalTime += float( clock () - assembleStar_BeginTime ) /  CLOCKS_PER_SEC;

      // FIND CORRELATION
      correlateStar_BeginTime = clock();
      corrCoeff = EvalCorrelationCoefficient(resVec,totalStarFaces,facesID,facesCoeffs);
      correlateStar_TotalTime += float( clock () - correlateStar_BeginTime ) /  CLOCKS_PER_SEC;

      // Store Correlation coefficient in Expansion
      if(!isBC){
        expansion->vortexCoeff[componentCount] += corrCoeff;
      }else{
        bcExpansion->vortexCoeff[componentCount] += corrCoeff;
      }

      // UPDATE RESIDUAL
      updateStar_BeginTime = clock();
      UpdateResidualAndFilter(corrCoeff,totalStarFaces,facesID,facesCoeffs,resVec,filteredVec,resNorm);
      updateStar_TotalTime += float( clock () - updateStar_BeginTime ) /  CLOCKS_PER_SEC;
    }

    // Eval Two-Norm of the Coefficient Vector
    if(!isBC){
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
      if(oldResNorm<kMathZero){
        converged = true;
      }else{
        converged = (fabs((resNorm-oldResNorm)/(oldResNorm))<itTol);
      }
    }else{
      converged = false;
    }

    // Update Norm
    oldResNorm = resNorm;
    oldTwoNorm = twoNorm;
  }

  // Iterations Stopped By User
  // TODO: ALLOW THE USER TO STOP THE ITERATIONS!!!
  if(stopIterations){
    WriteSchMessage("ITERATIONS STOPPED BY USER.");
  }

  // PRINT TIME STATISTICS
  printf("--- TIME STATISTICS\n");
  printf("Residual Assembly Time: %f [s]\n",assembleRes_TotalTime);
  printf("Constant Pattern Correlation Time: %f [s]\n",constPattern_TotalTime);
  printf("Star Shape Assembly Time: %f [s]\n",assembleStar_TotalTime);
  printf("Star Correlation Time: %f [s]\n",correlateStar_TotalTime);
  printf("Residual Update Time: %f [s]\n",updateStar_TotalTime);

  // Final Residual
  printf("--- INFO\n");
  WriteSchMessage("Final Residual Norm: " + MRIUtils::FloatToStr(resNorm) + "\n");

  // Check If the Flux Is Locally Conservative
  maxDivergence = EvalMaxDivergence(filteredVec);

  // Write Divergence Message
  WriteSchMessage("Max Divergence: " + MRIUtils::FloatToStr(maxDivergence) + "\n");

  // Make Diffence between Coefficient Expansions
  if(isBC){
    for(int loopA=0;loopA<3;loopA++){
      expansion->constantFluxCoeff[loopA] -= bcExpansion->constantFluxCoeff[loopA];
    }
    for(int loopA=0;loopA<totalVortexes;loopA++){
      expansion->vortexCoeff[loopA] -= bcExpansion->vortexCoeff[loopA];
    }
  }

  // Recover Velocities from Face Fluxes
  RecoverCellVelocitiesRT0(isBC,filteredVec);

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

  // CLEAR
  facesID.clear();
  facesCoeffs.clear();

  // RESIZE
  facesID.resize(10);
  facesCoeffs.resize(10);

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
  UpdateVelocities();

  // DELETE ARRAY
  delete [] faceFluxVec;
}
