#include "mriScan.h"
#include "mriUtils.h"

// ===========
// CONSTRUCTOR
// ===========
MRIScan::MRIScan(double currentTime){
  // Assign Scan Time
  scanTime = currentTime;
  // Resent Pressure Gradient and Relative Pressure
  hasPressureGradient = false;
  hasRelativePressure = false;
  hasReynoldsStress = false;
}

// ================
// COPY CONSTRUCTOR
// ================
MRIScan::MRIScan(const MRIScan &copyScan){
  // Assign Scan Time
  scanTime = copyScan.scanTime;
  // Resent Pressure Gradient and Relative Pressure
  hasPressureGradient = copyScan.hasPressureGradient;
  hasRelativePressure = copyScan.hasRelativePressure;
  hasReynoldsStress = copyScan.hasReynoldsStress;
  // Copy cells totals
  for(int loopA=0;loopA<3;loopA++){
    domainSizeMin[loopA] = copyScan.domainSizeMin[loopA];
    domainSizeMax[loopA] = copyScan.domainSizeMax[loopA];
  }
  // Max Velocity Module
  maxVelModule = 0.0;
  // Total number of Cells
  totalCellPoints = copyScan.totalCellPoints;
  // Allocate cell points
  cellPoints.resize(totalCellPoints);
  // Initialize
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Copy position
    cellPoints[loopA].position[0] = copyScan.cellPoints[loopA].position[0];
    cellPoints[loopA].position[1] = copyScan.cellPoints[loopA].position[1];
    cellPoints[loopA].position[2] = copyScan.cellPoints[loopA].position[2];
    // Concentration
    cellPoints[loopA].concentration = 0.0;
    // Velocity
    cellPoints[loopA].velocity[0] = 0.0;
    cellPoints[loopA].velocity[1] = 0.0;
    cellPoints[loopA].velocity[2] = 0.0;
    // Filtered Velocity
    cellPoints[loopA].auxVector[0] = 0.0;
    cellPoints[loopA].auxVector[1] = 0.0;
    cellPoints[loopA].auxVector[2] = 0.0;
    // Pressure Gradient
    cellPoints[loopA].pressGrad[0] = 0.0;
    cellPoints[loopA].pressGrad[1] = 0.0;
    cellPoints[loopA].pressGrad[2] = 0.0;
    // Relative Pressure
    cellPoints[loopA].relPressure = 0.0;
  }
}

// =================================
// UPDATE VELOCITIES FROM AUX VECTOR
// =================================
void MRIScan::UpdateVelocities(){
  maxVelModule = 0.0;
  MRIReal currentNorm = 0.0;
  // Update Velocities
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Assign Filtered Vectors
    cellPoints[loopA].velocity[0] = cellPoints[loopA].auxVector[0];
    cellPoints[loopA].velocity[1] = cellPoints[loopA].auxVector[1];
    cellPoints[loopA].velocity[2] = cellPoints[loopA].auxVector[2];
    // Get New Norm
    currentNorm = MRIUtils::Do3DEucNorm(cellPoints[loopA].auxVector);
    if(currentNorm>maxVelModule){
      maxVelModule = currentNorm;
    }
  }
}

// ====================
// EXPORT NODES TO FILE
// ====================
void MRIScan::ExportNodesToFile(std::string fileName){
    // Open Output File
    FILE* outFile;
    outFile = fopen(fileName.c_str(),"w");
    // Write Header
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    fprintf(outFile,"%e;%e;%e\n",cellPoints[loopA].position[0],cellPoints[loopA].position[1],cellPoints[loopA].position[2]);
  }
    // Close Output file
    fclose(outFile);
}

// ===================
// GET DIFFERENCE NORM
// ===================
double MRIScan::GetDiffNorm(MRIScan* otherScan){
  double diffNorm = 0;
  double currDiff[3] = {0.0};
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    currDiff[0] = cellPoints[loopA].velocity[0] - otherScan->cellPoints[loopA].velocity[0];
    currDiff[1] = cellPoints[loopA].velocity[1] - otherScan->cellPoints[loopA].velocity[1];
    currDiff[2] = cellPoints[loopA].velocity[2] - otherScan->cellPoints[loopA].velocity[2];
    diffNorm += sqrt(currDiff[0]*currDiff[0] + currDiff[1]*currDiff[1] + currDiff[2]*currDiff[2]);
  }
  return diffNorm;
}

// ===================================
// EVAL NOISY PRESSURE GRADIENT POINTS
// ===================================
void MRIScan::EvalNoisyPressureGradientPoints(){
  // DECLARE
  std::vector<int> neighbours;
  double currPressGradX = 0.0;
  double currPressGradY = 0.0;
  double currPressGradZ = 0.0;
  double avPressGradX = 0.0;
  double avPressGradY = 0.0;
  double avPressGradZ = 0.0;
  int currCell = 0;
  double currentDistance = 0.0;
  double currentDistance1 = 0.0;
  double currentDistance2 = 0.0;
  double currentDistance3 = 0.0;
  double currentModulus = 0.0;
  WriteSchMessage(std::string("Seeking Noisy Pressure Gradient Locations...\n"));
  // LOOP ON ALL INNER CELLS
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    if(IsInnerCell(loopA)){
      // GET PRESSURE GRADIENT AT CURRENT POINT
      currPressGradX = cellPoints[loopA].pressGrad[0];
      currPressGradY = cellPoints[loopA].pressGrad[1];
      currPressGradZ = cellPoints[loopA].pressGrad[2];
      // GET NEIGHBOURS
      // Should I use all neighbors including the current cell!!!
      GetCartesianNeighbourCells(loopA,neighbours,true);
      // INIT AVERAGE PRESSURE GRADIENT
      avPressGradX = 0.0;
      avPressGradY = 0.0;
      avPressGradZ = 0.0;
      // GET THE VALUES ON NEIGHBOR CELLS
      for(int loopB=0;loopB<k3DNeighbors;loopB++){
        currCell = neighbours[loopB];
        // Get Quantity for Neighbor Cell
        avPressGradX += cellPoints[currCell].pressGrad[0];
        avPressGradY += cellPoints[currCell].pressGrad[1];
        avPressGradZ += cellPoints[currCell].pressGrad[2];
      }
      // DIVIDE BY THE NUMBER OF COMPONENTS
      avPressGradX = avPressGradX/(double)6.0;
      avPressGradY = avPressGradY/(double)6.0;
      avPressGradZ = avPressGradZ/(double)6.0;
      // SET VALUE
      currentDistance = sqrt((currPressGradX-avPressGradX)*(currPressGradX-avPressGradX)+
                             (currPressGradY-avPressGradY)*(currPressGradY-avPressGradY)+
                             (currPressGradZ-avPressGradZ)*(currPressGradZ-avPressGradZ));

      currentDistance1 = fabs(currPressGradX-avPressGradX);
      currentDistance2 = fabs(currPressGradY-avPressGradY);
      currentDistance3 = fabs(currPressGradZ-avPressGradZ);
      // GET MODULUS
      currentModulus = sqrt((currPressGradX)*(currPressGradX)+
                             (currPressGradY)*(currPressGradY)+
                             (currPressGradZ)*(currPressGradZ));

      // ASSIGN VALUE AS CONCENTRATION
      //if(fabs(currentModulus)>1500.0){
      if(fabs(currentModulus)>1.0e-7){
        // JET JULIEN
        cellPoints[loopA].auxVector[0] = (currentDistance/currentModulus);
        //cellPoints[loopA].filteredVel[0] = currentModulus;
        //cellPoints[loopA].filteredVel[0] = (currentDistance1/currentModulus);
        //cellPoints[loopA].filteredVel[1] = (currentDistance2/currentModulus);
        //cellPoints[loopA].filteredVel[2] = (currentDistance3/currentModulus);
        //cellPoints[loopA].filteredVel[0] = currentModulus;
      }else{
        cellPoints[loopA].auxVector[0] = 0.0;
        //cellPoints[loopA].filteredVel[0] = 0.0;
        //cellPoints[loopA].filteredVel[1] = 0.0;
        //cellPoints[loopA].filteredVel[2] = 0.0;
      }
    }
  }
}

// =========================
// COMPUTE QUANTITY GRADIENT
// =========================
void MRIScan::ComputeQuantityGradient(int qtyID){
  hasPressureGradient = true;
  double* gradient = new double[3];
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    EvalSpaceGradient(loopA,qtyID,gradient);
    cellPoints[loopA].setQuantity(kQtyPressGradientX,gradient[0]);
    cellPoints[loopA].setQuantity(kQtyPressGradientY,gradient[1]);
    cellPoints[loopA].setQuantity(kQtyPressGradientZ,gradient[2]);
  }
}

// ==========================
// SAVE QUANTITIES TO OUTPUTS
// ==========================
void MRIScan::saveVelocity(){
  MRIOutput out1("Original_Velocity",3);
  double totalVel = 0.0;
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    out1.values.push_back(cellPoints[loopA].velocity[0]);
    out1.values.push_back(cellPoints[loopA].velocity[1]);
    out1.values.push_back(cellPoints[loopA].velocity[2]);
    totalVel += (cellPoints[loopA].velocity[0]*cellPoints[loopA].velocity[0] +
                 cellPoints[loopA].velocity[1]*cellPoints[loopA].velocity[1] +
                 cellPoints[loopA].velocity[2]*cellPoints[loopA].velocity[2]);
  }
  // Add output quantities
  outputs.push_back(out1);
  printf("Velocity Energy: %e\n",totalVel);
}

// ===============
// MESSAGE PASSING
// ===============
void MRIScan::DistributeScanData(MRICommunicator* comm){
  throw MRIException("ERROR: MRIScan::DistributeScanData Not Implemented.");
}
