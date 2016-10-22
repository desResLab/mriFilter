# include "mriScan.h"

void writeVectorToFile(string outFile, int size, double* vec){
  // Open Output File
  FILE* f;
  f = fopen(outFile.c_str(),"w");
  for(int loopA=0;loopA<size;loopA++){
    fprintf(f,"%f\n",vec[loopA]);
  }
  // Close Output file
  fclose(f);
}

// ===================
// PRINT FASE ID INDEX
// ===================
void printFaceIDIndexes(std::string fileName, int totalStarFaces, std::vector<int> facesID, std::vector<double> facesCoeffs){
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
void printResidualVector(std::string fileName, int totalFaces, double* resVec){
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

// =======================
// ASSEMBLE CONSTANT SHAPE
// =======================
void MRIScan::assembleConstantPattern(int currentDim, int& totalConstantFaces,
                                      MRIIntVec& facesID, MRIDoubleVec& facesCoeffs){

  // Clear Vectors
  totalConstantFaces = 0;
  facesID.clear();
  facesCoeffs.clear();
  MRIIntVec orientation;

  // Loop Over Faces
  int totNegFaces = 0;
  for(size_t loopA=0;loopA<topology->faceConnections.size();loopA++){
    if(fabs(topology->faceNormal[loopA][currentDim])>kMathZero){
      facesID.push_back(loopA);
      if(topology->faceNormal[loopA][currentDim]>kMathZero){
        orientation.push_back(1);
      }else{
        orientation.push_back(-1);
        totNegFaces++;
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

// =======================
// ASSEMBLE CONSTANT SHAPE
// =======================
void MRIScan::assembleConstantPatternMPI(int currentDim, int& totalConstantFacesOnProc,
                                         MRIIntVec& facesIDOnProc, MRIDoubleVec& facesCoeffsOnProc,
                                         int minFaceOnProc, int maxFaceOnProc, MRICommunicator* comm){

  // Clear Vectors
  totalConstantFacesOnProc = 0;
  facesIDOnProc.clear();
  facesCoeffsOnProc.clear();
  MRIDoubleVec orientation;

  // Get Total Number of Faces in this Direction
  int totFacesThisDir = 0;
  switch(currentDim){
    case 0:
      totFacesThisDir = topology->cellTotals[1] * topology->cellTotals[2] * (topology->cellTotals[0] + 1);
      break;
    case 1:
      totFacesThisDir = topology->cellTotals[0] * topology->cellTotals[2] * (topology->cellTotals[1] + 1);
      break;
    case 2:
      totFacesThisDir = topology->cellTotals[0] * topology->cellTotals[1] * (topology->cellTotals[2] + 1);
      break;
  }

  // Decide Orientations
  int totNegative = 0;
  for(int loopA=minFaceOnProc;loopA<maxFaceOnProc;loopA++){
    if(fabs(topology->faceNormal[loopA][currentDim])>kMathZero){
      facesIDOnProc.push_back(loopA);
      if(topology->faceNormal[loopA][currentDim]>kMathZero){
        orientation.push_back(1.0);
      }else{
        orientation.push_back(-1.0);
        totNegative++;
      }
    }
  }

  // Fill Coefficients
  for(size_t loopA=0;loopA<facesIDOnProc.size();loopA++){
    facesCoeffsOnProc.push_back((orientation[loopA]/sqrt((double)totFacesThisDir)));
  }

  // Update Counter
  totalConstantFacesOnProc = facesIDOnProc.size();
}


// ====================
// ASSEMBLE STAR MATRIX
// ====================
void MRIScan::assembleStarMatrix(int &totalFaces, int &totalBasis, MRIDoubleMat& starMatrix){
  // Init Rows and Columns
  totalFaces = topology->getTotalFaces();
  totalBasis = getTotalBasisNumber();
  // Init Count
  int totalStarFaces = 0;
  std::vector<int> facesID;
  std::vector<double> facesCoeffs;
  int currentRow = 0;
  int currentCount = 0;
  // Get total Number of Vortexes
  int totalVortices = evalTotalVortex();
  // Allocate and Initialize
  starMatrix.resize(totalFaces);
  for(int loopA=0;loopA<totalFaces;loopA++){
    starMatrix[loopA].resize(totalBasis); 
  }
  // Loop On The Three Directions
	int totalSlices = 0;
	int totalStars = 0;
  for(int loopB=0;loopB<totalVortices;loopB++){
    // Find Star Shape
    assembleStarShape(loopB,totalStarFaces,facesID,facesCoeffs);
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
int MRIScan::getTotalBasisNumber(){
  int totBasis = 0;
  int totalSlices = 0;
  int totalStars = 0;
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    switch(loopA){
      case 0:
        // YZ Planes
        totalSlices = topology->cellTotals[0];
        totalStars = (topology->cellTotals[1] + 1)*(topology->cellTotals[2] + 1);
        break;
      case 1:
        // XZ Planes
        totalSlices = topology->cellTotals[1];
        totalStars = (topology->cellTotals[0] + 1)*(topology->cellTotals[2] + 1);
        break;
      case 2:
        // XY Planes
        totalSlices = topology->cellTotals[2];
        totalStars = (topology->cellTotals[0] + 1)*(topology->cellTotals[1] + 1);
        break;
    }
    totBasis = totBasis + totalSlices * totalStars;
  }
	return totBasis;
}

// =================
// CHECK PERMUTATION
// =================
bool checkPermutation(int size, int* perm){
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
void MRIScan::expandStarShape(int totalStarFaces, int* facesID, double* facesCoeffs, double* &fullStarVector){
  // Get Total Number Of Faces
  int totalFaces = topology->cellTotals[0] * topology->cellTotals[1] * (topology->cellTotals[2] + 1)+
                   topology->cellTotals[1] * topology->cellTotals[2] * (topology->cellTotals[0] + 1)+
                   topology->cellTotals[2] * topology->cellTotals[0] * (topology->cellTotals[1] + 1);
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
double MRIScan::evalMaxDivergence(const MRIDoubleVec& filteredVec){
  double maxDivergence = 0.0;
  double currentDiv = 0.0;
  int currFace = 0;
  MRIDoubleVec extNormal(3,0.0);
  double normalSign = 0.0;
  // Loop on cells
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    currentDiv = 0.0;
    // Loop on faces
    for(size_t loopB=0;loopB<topology->cellFaces[loopA].size();loopB++){
      // Get Current Face
      currFace = topology->cellFaces[loopA][loopB];
      // Get External Normal
      topology->getExternalFaceNormal(loopA,loopB,extNormal);
      // Get Sign
      normalSign = 0.0;
      for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
        normalSign += extNormal[loopC] * topology->faceNormal[currFace][loopC];
      }
      MRIUtils::round(normalSign);
      // Get Sign
      currentDiv += filteredVec[topology->cellFaces[loopA][loopB]] * normalSign;
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
void MRIScan::evalCellDivergences(const MRIDoubleVec& faceVec, MRIDoubleVec& cellDivs){
  double currentDiv = 0.0;
  int currFace = 0;
  MRIDoubleVec extNormal(3,0.0);
  double normalSign = 0.0;

  // Loop on cells
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    currentDiv = 0.0;
    // Loop on faces
    for(size_t loopB=0;loopB<topology->cellFaces[loopA].size();loopB++){
      // Get Current Face
      currFace = topology->cellFaces[loopA][loopB];
      // Get External Normal
      topology->getExternalFaceNormal(loopA,loopB,extNormal);
      // Get Sign
      normalSign = 0.0;
      for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
        normalSign += extNormal[loopC] * topology->faceNormal[currFace][loopC];
      }
      MRIUtils::round(normalSign);
      // Get Sign
      currentDiv += faceVec[topology->cellFaces[loopA][loopB]] * normalSign;
    }
    // Store Divergence
    cellDivs.push_back(currentDiv);
  }
}

// ===================================
// RECOVER VELOCITIES FROM FACE FLUXES
// ===================================
void MRIScan::recoverCellVelocitiesRT0(bool useBCFilter, MRIDoubleVec& filteredVec){
  // Variables
  int locFace1 = 0;
  int locFace2 = 0;
  int firstFaceID = 0;
  int secondFaceID = 0;
  double firstAvVelocity = 0.0;
  double secondAvVelocity = 0.0;
  // Loop To Assemble Residual Vector
  for(int loopA=0;loopA<topology->totalCells;loopA++){
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
      firstFaceID = topology->cellFaces[loopA][locFace1];
      secondFaceID = topology->cellFaces[loopA][locFace2];
      // Eval Average Velocity
      firstAvVelocity = (filteredVec[firstFaceID] * topology->faceNormal[firstFaceID][loopB]/topology->faceArea[firstFaceID]);
      secondAvVelocity = (filteredVec[secondFaceID]* topology->faceNormal[secondFaceID][loopB]/topology->faceArea[secondFaceID]);
      // Set Correction
      if(useBCFilter){
        cells[loopA].auxVector[loopB] = cells[loopA].velocity[loopB]-0.5*(firstAvVelocity + secondAvVelocity);
      }else{
        cells[loopA].auxVector[loopB] = 0.5*(firstAvVelocity + secondAvVelocity);
      }
    }
  }
}

// =======================================================
// GET DIMENSION, SLICE AND STAR NUMBER FROM VORTEX NUMBER
// =======================================================
void MRIScan::getDimensionSliceStarFromVortex(int vortexNumber,int &dimNumber,int &sliceNumber,int &starNumber){
  // Declare
  int totalSlices[3] = {0};
  int totalStars[3] = {0};
  int vortexOffset = 0;
  
  // Get Dimension,Slice and Star from the cardinality
  // X Direction
  totalSlices[0] = topology->cellTotals[0];
  totalStars[0] = (topology->cellTotals[1]+1)*(topology->cellTotals[2]+1);
  // Y Direction
  totalSlices[1] = topology->cellTotals[1];
  totalStars[1] = (topology->cellTotals[0]+1)*(topology->cellTotals[2]+1);
  // Z Direction
  totalSlices[2] = topology->cellTotals[2];
  totalStars[2] = (topology->cellTotals[0]+1)*(topology->cellTotals[1]+1);

  // Get the total Vortexes for every dimension
  double totalVortexDim1 = totalSlices[0] * totalStars[0];
  double totalVortexDim2 = totalVortexDim1 + totalSlices[1] * totalStars[1];
  double totalVortexDim3 = totalVortexDim2 + totalSlices[2] * totalStars[2];

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
void MRIScan::assembleStarShape(int vortexNumber, int &totalFaces, std::vector<int> &facesID, std::vector<double> &facesCoeffs){
  // Declare
  int currFace = 0;
  double currCoeff = 0.0;

  // Clear Array and Reset Counters
  totalFaces = 0;
  //facesID.reserve(10);
  //facesCoeffs.reserve(10);

  // Loop Over All Faces
  double norm = 0.0;
  for(size_t loopB=0;loopB<topology->edgeFaces[vortexNumber].size();loopB++){
    currFace = abs(topology->edgeFaces[vortexNumber][loopB])-1;
    currCoeff = (double)topology->edgeFaces[vortexNumber][loopB]/(double)fabs(topology->edgeFaces[vortexNumber][loopB]);

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

// ===========================================
// UPDATE THE RESIDUAL AND RESIDUAL NORM - MPI
// ===========================================
void updateResidualAndFilterMPI(MRICommunicator* comm,
                                int globTotFaces, double corrCoeff, 
                                int totalFaces, 
                                const MRIIntVec& facesID, 
                                const MRIDoubleVec& facesCoeffs,                                
                                const MRIIntVec& minFaceGlob,
                                const MRIIntVec& maxFaceGlob,
                                double &resNorm,
                                MRIDoubleVec& resVec, 
                                MRIDoubleVec& filteredVels){

  double normSqrIncr = 0.0;
  double normSqrTot = 0.0;
  int currentFaceID = 0;
  double currentFaceCoeff = 0.0;
  for(int loopA=0;loopA<totalFaces;loopA++){
    currentFaceID = facesID[loopA];
    currentFaceCoeff = facesCoeffs[loopA];
    // Update Residual Norm
    normSqrIncr += (resVec[currentFaceID]-corrCoeff*currentFaceCoeff)*(resVec[currentFaceID]-corrCoeff*currentFaceCoeff) - (resVec[currentFaceID]*resVec[currentFaceID]);
    // UPDATE RESIDUAL AND RECONSTRUCTED FIELD
    // Update Residual
    resVec[currentFaceID] = resVec[currentFaceID] - corrCoeff * currentFaceCoeff;
    // Update Filtered Vector
    filteredVels[currentFaceID] = filteredVels[currentFaceID] + corrCoeff * currentFaceCoeff;
  }

  //printf("[%d] NORMINCR BEFORE All Reduce: %f\n",comm->currProc,normSqrIncr);

  // Reduce Norm Increment
  MPI_Allreduce(&normSqrIncr,&normSqrTot,1,MPI_DOUBLE,MPI_SUM,comm->mpiComm);

  //printf("[%d] NORMINCR AFTER All Reduce: %f\n",comm->currProc,normSqrTot);

  // Update Residual Norm
  resNorm = (resNorm*resNorm) + normSqrTot;
  resNorm = sqrt(fabs(resNorm));

  // Communicate the residual and Filtered Vector
  int size = maxFaceGlob[comm->currProc] - minFaceGlob[comm->currProc];
  int* recvcounts = new int[comm->totProc];
  int* displs = new int[comm->totProc];
  double* storeVec = new double[size];
  for(int loopA=0;loopA<comm->totProc;loopA++){
    displs[loopA] = minFaceGlob[loopA];
    recvcounts[loopA] = maxFaceGlob[loopA] - minFaceGlob[loopA];
  }

  // Residual Vector
  for(int loopA=minFaceGlob[comm->currProc];loopA<maxFaceGlob[comm->currProc];loopA++){
    storeVec[loopA-minFaceGlob[comm->currProc]] = resVec[loopA];
  }
  MPI_Allgatherv(&storeVec[0],size,MPI_DOUBLE,&resVec[0],recvcounts,displs,MPI_DOUBLE,comm->mpiComm);
  // Filter Vector
  for(int loopA=minFaceGlob[comm->currProc];loopA<maxFaceGlob[comm->currProc];loopA++){
    storeVec[loopA-minFaceGlob[comm->currProc]] = filteredVels[loopA];
  }
  MPI_Allgatherv(&storeVec[0],size,MPI_DOUBLE,&filteredVels[0],recvcounts,displs,MPI_DOUBLE,comm->mpiComm);
  // Free Memory
  delete [] recvcounts;
  delete [] displs;
  delete [] storeVec;
}

// =====================================
// UPDATE THE RESIDUAL AND RESIDUAL NORM
// =====================================
void updateResidualAndFilter(double corrCoeff, 
                             int totalFaces, 
                             const MRIIntVec& facesID, 
                             const MRIDoubleVec& facesCoeffs,
                             double &resNorm,
                             MRIDoubleVec& resVec, 
                             MRIDoubleVec& filteredVels){

  double normSqrIncr = 0.0;
  int currentFaceID = 0;
  double currentFaceCoeff = 0.0;
  for(int loopA=0;loopA<totalFaces;loopA++){
    currentFaceID = facesID[loopA];
    currentFaceCoeff = facesCoeffs[loopA];
    // Update Residual Norm
    normSqrIncr += (resVec[currentFaceID]-corrCoeff*currentFaceCoeff)*(resVec[currentFaceID]-corrCoeff*currentFaceCoeff) - (resVec[currentFaceID]*resVec[currentFaceID]);
    // UPDATE RESIDUAL AND RECONSTRUCTED FIELD
    // Update Residual
    resVec[currentFaceID] = resVec[currentFaceID] - corrCoeff * currentFaceCoeff;
    // Update Filtered Vector
    filteredVels[currentFaceID] = filteredVels[currentFaceID] + corrCoeff * currentFaceCoeff;
  }

  // Update Residual Norm
  resNorm = (resNorm*resNorm) + normSqrIncr;
  resNorm = sqrt(fabs(resNorm));
}



// EVAL CORRELATION COEFFICIENT
double evalCorrelationCoefficient(const MRIDoubleVec& resVec, int totalFaces, const MRIIntVec& facesID, const MRIDoubleVec& facesCoeffs){
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
std::string getFaceLocationString(int faceLocation){
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
void MRIScan::assembleResidualVector(bool useBCFilter, MRIThresholdCriteria* thresholdCriteria,
                                     int& totalFaces, MRIDoubleVec& resVec, MRIDoubleVec& filteredVec, double& resNorm){
  bool   continueToProcess = false;
  double currentValue = 0.0;
  double faceComponent = 0.0;
  int    currentFace = 0;
  double currFaceArea = 0.0;
  bool   checkPassed = false;

  // Get Total Number Of Faces
  totalFaces = topology->faceConnections.size();

  // Allocate
  resVec.resize(totalFaces);
  filteredVec.resize(totalFaces);
  MRIIntVec resID(totalFaces);

  // Initialize
  for(int loopA=0;loopA<totalFaces;loopA++){
    resVec[loopA] = 0.0;
    resID[loopA] = 0;
    filteredVec[loopA] = 0.0;
  }

  // Loop To Assemble Residual Vector
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Check for BC
    if(useBCFilter){
      currentValue = cells[loopA].getQuantity(thresholdCriteria->thresholdQty);
      continueToProcess = thresholdCriteria->meetsCriteria(currentValue);
    }else{
      continueToProcess = true;
    }
    if(continueToProcess){
      // Loop On Faces
      for(int loopB=0;loopB<k3DNeighbors;loopB++){
        // Get Current Face
        currentFace = topology->cellFaces[loopA][loopB];
        // Get Face Area
        currFaceArea = topology->faceArea[currentFace];
        // Get Normal Veclocity
        faceComponent = 0.0;
        for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
          faceComponent += cells[loopA].velocity[loopC] * topology->faceNormal[currentFace][loopC];
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
      std::string currentMsgs = "Internal: Wrong Face Connectivity, Face: " + MRIUtils::intToStr(loopA)+ "; Connectivity: " + MRIUtils::intToStr(resID[loopA])+".";
      throw new MRIException(currentMsgs.c_str());
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
}


// ================
// FORM VORTEX LIST
// ================
void MRIScan::formVortexList(MRICommunicator* comm,
                             int totVortex,
                             const MRIIntVec& minFace,
                             const MRIIntVec& maxFace,
                             MRIIntVec& innerVortexList,
                             MRIIntVec& boundaryVortexList){
  bool foundFace = true;
  bool isInternal = true;
  int currFace = 0;
  int countFace = 0;

  // Form Inner List for Current Processor
  for(int loopA=0;loopA<totVortex;loopA++){
    foundFace = true;
    for(size_t loopB=0;loopB<topology->edgeFaces[loopA].size();loopB++){
      currFace = abs(topology->edgeFaces[loopA][loopB])-1;
      foundFace = (foundFace && ((currFace>=minFace[comm->currProc])&&(currFace<maxFace[comm->currProc])));
    }
    if(foundFace){
      innerVortexList.push_back(loopA);
      countFace++;
    }
  }

  // Form Boundary List
  for(int loopA=0;loopA<totVortex;loopA++){
    foundFace = true;
    isInternal = false;
    for(int loopC=0;loopC<comm->totProc;loopC++){
      foundFace = true;
      for(size_t loopB=0;loopB<topology->edgeFaces[loopA].size();loopB++){
        currFace = abs(topology->edgeFaces[loopA][loopB])-1;
        foundFace = (foundFace && ((currFace >= minFace[loopC])&&(currFace < maxFace[loopC])));
      }
      isInternal = (isInternal || foundFace);
    }
    if(!isInternal){
      boundaryVortexList.push_back(loopA);
    }
  }
}

// =========================================
// COMMUNICATE RESIDUAL AND FILTERED VECTORS
// =========================================
void communicateResFilt2(MRICommunicator* comm,
                         int totalFaces,
                         const MRIIntVec& minFaceGlob,
                         const MRIIntVec& maxFaceGlob,
                         MRIDoubleVec& resVec,
                         MRIDoubleVec& filteredVels,
                         double& resNorm){
  
  int size = maxFaceGlob[comm->currProc] - minFaceGlob[comm->currProc];
  MRIIntVec recvcounts(comm->totProc);
  MRIIntVec displs(comm->totProc);
  MRIDoubleVec storeVec(size);
  
  for(int loopA=0;loopA<comm->totProc;loopA++){
    displs[loopA] = minFaceGlob[loopA];
    recvcounts[loopA] = maxFaceGlob[loopA] - minFaceGlob[loopA];
  }
  
  // Residual Vector
  for(int loopA=minFaceGlob[comm->currProc];loopA<maxFaceGlob[comm->currProc];loopA++){
    storeVec[loopA-minFaceGlob[comm->currProc]] = resVec[loopA];
  }
  MPI_Allgatherv(&storeVec[0],size,MPI_DOUBLE,&resVec[0],&recvcounts[0],&displs[0],MPI_DOUBLE,comm->mpiComm);
  
  // Filter Vector
  for(int loopA=minFaceGlob[comm->currProc];loopA<maxFaceGlob[comm->currProc];loopA++){
    storeVec[loopA-minFaceGlob[comm->currProc]] = filteredVels[loopA];
  }
  MPI_Allgatherv(&storeVec[0],size,MPI_DOUBLE,&filteredVels[0],&recvcounts[0],&displs[0],MPI_DOUBLE,comm->mpiComm);
  // Update Norm
  resNorm = 0.0;
  for(int loopA=0;loopA<totalFaces;loopA++){
    resNorm += resVec[loopA] * resVec[loopA];
  }
  resNorm = sqrt(resNorm);
}

// =========================================
// SYNC TEMPORARY EXPANSION AMONG PROCESSORS
// =========================================
void syncExpansion(int totalVortexes, MRIDoubleVec& currentExp, MRICommunicator* comm){
  MRIDoubleVec target(totalVortexes);
  MPI_Reduce(&currentExp[0],&target[0],totalVortexes,MPI_DOUBLE,MPI_SUM,0,comm->mpiComm);
  if(comm->currProc == 0){
    for(int loopA=0;loopA<totalVortexes;loopA++){
      currentExp[loopA] = target[loopA];
    }
  }
}

// =================
// PHYSICS FILTERING
// =================
void MRIScan::applySMPFilter(MRICommunicator* comm, bool isBC, 
                             MRIThresholdCriteria* thresholdCriteria,
                             double itTol,
                             int maxIt,
                             bool useConstantPatterns){

  // INITIALIZATION
  int totalFaces = topology->faceConnections.size();
  MRIDoubleVec resVec;
  MRIDoubleVec filteredVec;
  MRIIntVec facesID;
  MRIDoubleVec facesCoeffs;
  MRIIntVec innerVortexList;
  MRIIntVec boundaryVortexList;
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
  double redcorrCoeff = 0.0;

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

  // Global Array for Face Storage
  MRIIntVec minFaceGlob(comm->totProc);
  MRIIntVec maxFaceGlob(comm->totProc);
  int minFaceOnProc = 0.0;
  int maxFaceOnProc = 0.0;
  for(int loopA=0;loopA<comm->totProc;loopA++){
    // Determine the Minimum and Maximum Face number of current processor
    minFaceOnProc = int(loopA * (totalFaces - 1)/(comm->totProc));
    maxFaceOnProc = int((loopA + 1) * (totalFaces - 1)/(comm->totProc));
    if(loopA == (comm->totProc-1)){
      maxFaceOnProc = totalFaces;
    }
    // Store Face Limits
    minFaceGlob[loopA] = minFaceOnProc;
    maxFaceGlob[loopA] = maxFaceOnProc;
  }

  // Determine the Minimum and Maximum Face number of current processor
  minFaceOnProc = minFaceGlob[comm->currProc];
  maxFaceOnProc = maxFaceGlob[comm->currProc];

  // All processes are waiting for the root to read the files
  mpiError = MPI_Barrier(comm->mpiComm);
  MRIUtils::checkMpiError(mpiError);

  printf("[%d/%d] MinFace: %d, MaxFace: %d\n",comm->currProc,comm->totProc,minFaceOnProc,maxFaceOnProc);

  // Assemble Face Flux Vectors
  assembleRes_BeginTime = clock();
  assembleResidualVector(isBC,thresholdCriteria,totalFaces,resVec,filteredVec,resNorm);
  assembleRes_TotalTime += float( clock () - assembleRes_BeginTime ) /  CLOCKS_PER_SEC;

  // Initial Residual
  if(comm->currProc == 0){
    writeSchMessage("\n");
    if (isBC){
      writeSchMessage("FILTER ALGORITHM - BC - Step: "+MRIUtils::floatToStr(scanTime)+" ---------------------------\n");
    }else{
      writeSchMessage("FILTER ALGORITHM - FULL - Step "+MRIUtils::floatToStr(scanTime)+" ---------------------------\n");
    }
  }

  // START CLOCK
  const clock_t begin_time = clock();
  
  if(comm->currProc == 0){
    writeSchMessage("Initial Residual Norm: "+MRIUtils::floatToStr(resNorm)+"\n");
  }

  // Initialize Expansion
  MRIExpansion* bcExpansion = NULL;
  int totalVortexes = evalTotalVortex();

  // Allocate Temporary Expansion coefficient place holder
  MRIDoubleVec currentExp(totalVortexes);

  // Form List of Vortexes Including Faces for MPI
  if(comm->totProc > 1){
    formVortexList(comm,totalVortexes,minFaceGlob,maxFaceGlob,innerVortexList,boundaryVortexList);
  }else{
    for(int loopA=0;loopA<totalVortexes;loopA++){
      innerVortexList.push_back(loopA);
    }
  }

  // Processor 0 has expansion Coefficients
  if(comm->currProc == 0){
    if(!isBC){
      expansion = new MRIExpansion(totalVortexes);
    }else{
      bcExpansion = new MRIExpansion(totalVortexes);
    }
  }

  // Apply MP Filter
  bool converged = false;
  int itCount = 0;
  double oldResNorm = resNorm;
  double oldTwoNorm = twoNorm;
  // Initialize Component Count
  int componentCount = 0;

  // Start Filter Loop
  while((!converged)&&(itCount<maxIt)){

    // Update Iteration Count
    itCount++;

    constPattern_BeginTime = clock();

    // LOOP ON THE THREE DIRECTIONS
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      // Correlate Constant Patterns
      if(useConstantPatterns){

        // Find Constant Patterns
        if(comm->totProc > 1){
          assembleConstantPatternMPI(loopB,totalStarFaces,facesID,facesCoeffs,minFaceOnProc,maxFaceOnProc,comm);
        }else{
          assembleConstantPattern(loopB,totalStarFaces,facesID,facesCoeffs);
        }

        // Find Correlation
        corrCoeff = evalCorrelationCoefficient(resVec,totalStarFaces,facesID,facesCoeffs);

        // Store Correlation coefficient in Expansion
        if(comm->totProc > 1){
          MPI_Allreduce(&corrCoeff, &redcorrCoeff, 1, MPI_DOUBLE, MPI_SUM, comm->mpiComm);
          corrCoeff = redcorrCoeff;
        }

        // Store Expansion Coefficients on Master Processor
        if(comm->currProc == 0){
          if(!isBC){
            expansion->constantFluxCoeff[loopB] += corrCoeff;
          }else{
            bcExpansion->constantFluxCoeff[loopB] += corrCoeff;
          }
        }

        // Update Residual
        if(comm->totProc > 1){
          updateResidualAndFilterMPI(comm,totalFaces,corrCoeff,
                                     totalStarFaces,facesID,facesCoeffs,                                
                                     minFaceGlob,maxFaceGlob,
                                     resNorm,resVec,filteredVec);
        }else{
          updateResidualAndFilter(corrCoeff,totalStarFaces,facesID,facesCoeffs,resNorm,resVec,filteredVec);
        }
      }
    }

    //string fileName("faceEdge_proc_" + to_string(comm->currProc) + ".dat");
    //printfIntMatToFile(fileName,faceEdges);

    constPattern_TotalTime = float( clock () - constPattern_BeginTime ) /  CLOCKS_PER_SEC;

    // Reserve these two
    facesID.resize(10);
    facesCoeffs.resize(10);

    // LOOP ON VORTEXES
    componentCount = -1;

    int currVortex = 0;
    // Clean expansion
    for(int loopB=0;loopB<totalVortexes;loopB++){
      currentExp[loopB] = 0.0;
    }
    for(int loopB=0;loopB<innerVortexList.size();loopB++){
      // Increment the current component
      componentCount++;

      // Get Current Vortex
      currVortex = innerVortexList[loopB];

      // ASSEMBLE STAR
      assembleStar_BeginTime = clock();
      assembleStarShape(currVortex,totalStarFaces,facesID,facesCoeffs);
      // Communicate
      assembleStar_TotalTime += float( clock () - assembleStar_BeginTime ) /  CLOCKS_PER_SEC;

      // FIND CORRELATION
      correlateStar_BeginTime = clock();
      corrCoeff = evalCorrelationCoefficient(resVec,totalStarFaces,facesID,facesCoeffs);
      correlateStar_TotalTime += float( clock () - correlateStar_BeginTime ) /  CLOCKS_PER_SEC;

      // Store Correlation coefficient in Expansion
      currentExp[currVortex] += corrCoeff;

      // UPDATE RESIDUAL
      updateStar_BeginTime = clock();
      updateResidualAndFilter(corrCoeff,totalStarFaces,facesID,facesCoeffs,resNorm,resVec,filteredVec);
      updateStar_TotalTime += float( clock () - updateStar_BeginTime ) /  CLOCKS_PER_SEC;
    }

    // Sync Expansion for main and boundary filter
    if(comm->totProc > 1){
      syncExpansion(totalVortexes,currentExp,comm);
    }
    // Add to stored expansion
    if(comm->currProc == 0){
      for(int loopB=0;loopB<totalVortexes;loopB++){
        if(!isBC){
          expansion->vortexCoeff[loopB] += currentExp[loopB];
        }else{
          bcExpansion->vortexCoeff[loopB] += currentExp[loopB];
        }
      }
    }

    // IF MPI then Communicate Residual Vector
    if(comm->totProc > 1){
      // Communicate Residual, FilteredVels and Update Norm
      communicateResFilt2(comm,totalFaces,
                          minFaceGlob,maxFaceGlob,
                          resVec,filteredVec,resNorm);
    }

    for(int loopB=0;loopB<boundaryVortexList.size();loopB++){
      // Increment the current component
      componentCount++;

      // Get Current Vortex
      currVortex = boundaryVortexList[loopB];

      // ASSEMBLE STAR
      assembleStar_BeginTime = clock();
      assembleStarShape(currVortex,totalStarFaces,facesID,facesCoeffs);
      // Communicate
      assembleStar_TotalTime += float( clock () - assembleStar_BeginTime ) /  CLOCKS_PER_SEC;

      // FIND CORRELATION
      correlateStar_BeginTime = clock();
      corrCoeff = evalCorrelationCoefficient(resVec,totalStarFaces,facesID,facesCoeffs);
      correlateStar_TotalTime += float( clock () - correlateStar_BeginTime ) /  CLOCKS_PER_SEC;

      // Store Correlation coefficient in Expansion
      if(comm->currProc == 0){
        if(!isBC){
          expansion->vortexCoeff[currVortex] += corrCoeff;
        }else{
          bcExpansion->vortexCoeff[currVortex] += corrCoeff;
        }
      }

      // UPDATE RESIDUAL
      updateStar_BeginTime = clock();
      updateResidualAndFilter(corrCoeff,totalStarFaces,facesID,facesCoeffs,resNorm,resVec,filteredVec);
      updateStar_TotalTime += float( clock () - updateStar_BeginTime ) /  CLOCKS_PER_SEC;
    }

    //printf("[%d] RESIDUAL AFTER STARS: %f\n",comm->currProc,resNorm);

    // Eval Two-Norm of the Coefficient Vector
    if(comm->currProc == 0){
      if(!isBC){
        twoNorm = expansion->get2Norm(false);
      }else{
        twoNorm = bcExpansion->get2Norm(false);
      }
    }
    // Broadcase Two Norm
    MPI_Bcast(&twoNorm,1,MPI_DOUBLE,0,comm->mpiComm);

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
    if(comm->currProc == 0){
      writeSchMessage("[" + MRIUtils::intToStr(comm->currProc) + "] It: " + MRIUtils::intToStr(itCount) + "; ABS Res: "+MRIUtils::floatToStr(resNorm)+"; Rel: " + MRIUtils::floatToStr(relResNorm) +
                      "; Coeff 2-Norm: "+MRIUtils::floatToStr(twoNorm)+"; Rel 2-Norm: " + MRIUtils::floatToStr(relTwoNorm)+"\n");
    }

    // Check Convergence
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

  // Check If the Flux Is Locally Conservative
  maxDivergence = evalMaxDivergence(filteredVec);

  // Make Diffence between Coefficient Expansions
  if(comm->currProc == 0){
    if(isBC){
      for(int loopA=0;loopA<3;loopA++){
        expansion->constantFluxCoeff[loopA] -= bcExpansion->constantFluxCoeff[loopA];
      }
      for(int loopA=0;loopA<totalVortexes;loopA++){
        expansion->vortexCoeff[loopA] -= bcExpansion->vortexCoeff[loopA];
      }
    }
  }

  // Recover Velocities from Face Fluxes
  recoverCellVelocitiesRT0(isBC,filteredVec);

  // WRITE CPU TIME AND NUMBER OF ITERATIONS
  float totalCPUTime = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
  writeSchMessage("Total Iterations " + MRIUtils::intToStr(itCount) + "; Total CPU Time: " + MRIUtils::floatToStr(totalCPUTime) + "\n");

  // Eval Magniture and Angle Error
  recoverGlobalErrorEstimates(maxNormError,maxAngleError);

  // PRINT TIME STATISTICS
  if(comm->currProc == 0){
    printf("--- TIME STATISTICS\n");
    printf("Residual Assembly Time: %f [s]\n",assembleRes_TotalTime);
    printf("Constant Pattern Correlation Time: %f [s]\n",constPattern_TotalTime);
    printf("Star Shape Assembly Time: %f [s]\n",assembleStar_TotalTime);
    printf("Star Correlation Time: %f [s]\n",correlateStar_TotalTime);
    printf("Residual Update Time: %f [s]\n",updateStar_TotalTime);
    printf("\n");
    printf("--- INFO\n");
    writeSchMessage("Final Residual Norm: " + MRIUtils::floatToStr(resNorm) + "\n");
    // Write Divergence Message
    writeSchMessage("Max Divergence: " + MRIUtils::floatToStr(maxDivergence) + "\n");
    // Write Final Residual
    writeSchMessage("Average Magnitude Error: " + MRIUtils::floatToStr(maxNormError) + "\n");
    writeSchMessage("Average Angular Error: " + MRIUtils::floatToStr(maxAngleError) + "\n");
    // Close Filter Comments
    writeSchMessage("---\n");
  }

  // Update Velocities
  //UpdateVelocities(GlobalData);

  // Filtered Velocities Are Available
  //AreVelocityFiltered:=TRUE;
}

// ==================================
// REBUILD FROM EXPANSION COEFFICIENT
// ==================================
void MRIScan::rebuildFromExpansion(MRIExpansion* expansion, bool useConstantFlux){
  
  // RECONSTRUCTION
  int totalStarFaces = 0;
  int currFaceID = 0;
  double currFaceCoeff = 0.0;
  MRIIntVec facesID;
  MRIDoubleVec facesCoeffs;
  
  // Get total number of vortexes
  int totalVortex = evalTotalVortex();
  MRIDoubleVec faceFluxVec(totalVortex + 3);

  // INITIALIZE
  for(int loopA=0;loopA<totalVortex+3;loopA++){
    faceFluxVec[loopA] = 0.0;
  }

  // GLOBAL ATOMS
  if(useConstantFlux){
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      // Find Star Shape
      assembleConstantPattern(loopB/*Current Dimension*/,totalStarFaces,facesID,facesCoeffs);
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
    assembleStarShape(loopB,totalStarFaces,facesID,facesCoeffs);
    // Find Correlation
    for(int loopE=0;loopE<totalStarFaces;loopE++){
      currFaceID = facesID[loopE];
      currFaceCoeff = facesCoeffs[loopE];
      faceFluxVec[currFaceID] += currFaceCoeff*expansion->vortexCoeff[componentCount];
    }
  }

  // Recover Velocities from Face Fluxes
  recoverCellVelocitiesRT0(false,faceFluxVec);

  // Update Velocities
  updateVelocities();
}
