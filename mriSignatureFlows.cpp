#include <math.h>
#include "mriScan.h"
#include "mriConstants.h"
#include "schMessages.h"
#include "mriUtils.h"

// STAGNATION FLOW SOLUTION
void MRIScan::AssignStagnationFlowSignature(MRIDirection dir){
  double bConst = 1.0;
  double xCoord,yCoord,zCoord;
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Set to Zero
    cellPoints[loopA].concentration = 0.0;
    switch(dir){
      case kdirX:
        xCoord = cellPoints[loopA].position[0] - 0.0;
        yCoord = cellPoints[loopA].position[1] - domainSizeMin[1];
        zCoord = cellPoints[loopA].position[2] - 0.5*(domainSizeMin[2]+domainSizeMax[2]);
        cellPoints[loopA].velocity[0] = 0.0;
        cellPoints[loopA].velocity[1] = -bConst * yCoord;
        cellPoints[loopA].velocity[2] = bConst * zCoord;
        break;
      case kdirY:
        xCoord = cellPoints[loopA].position[0] - domainSizeMin[0];
        yCoord = cellPoints[loopA].position[1]- 0.0;
        zCoord = cellPoints[loopA].position[2]- 0.5 * (domainSizeMin[2] + domainSizeMax[2]);
        cellPoints[loopA].velocity[0] = bConst * xCoord;
        cellPoints[loopA].velocity[1] = 0.0;
        cellPoints[loopA].velocity[2] = -bConst * zCoord;
        break;
      case kdirZ:
        xCoord = cellPoints[loopA].position[0] - domainSizeMin[1];
        yCoord = cellPoints[loopA].position[1] - 0.5 * (domainSizeMin[1] + domainSizeMax[1]);
        zCoord = cellPoints[loopA].position[2] - 0.0;
        cellPoints[loopA].velocity[0] = bConst * xCoord;
        cellPoints[loopA].velocity[1] = -bConst * yCoord;
        cellPoints[loopA].velocity[2] = 0.0;
        break;
    }
  }
}

// EVAL TANGENT DIRECTION
void EvalTangentDirection(MRIDirection dir, double* radialVector, double* tangVector){
  double axialVec[3] = {0.0};
  switch(dir){
    case kdirX:
      axialVec[0] = 1.0;
      axialVec[1] = 0.0;
      axialVec[2] = 0.0;
      break;
    case kdirY:
      axialVec[0] = 0.0;
      axialVec[1] = 1.0;
      axialVec[2] = 0.0;
      break;
    case kdirZ:
      axialVec[0] = 0.0;
      axialVec[1] = 0.0;
      axialVec[2] = 1.0;
      break;    
  }
  // Perform External Product
  MRIUtils::Do3DExternalProduct(axialVec,radialVector,tangVector);
  // Normalize Result
  MRIUtils::Normalize3DVector(tangVector);
}

// ASSIGN CYLINDRICAL VORTEX FLOW 
void MRIScan::AssignCylindricalFlowSignature(MRIDirection dir){
  // Set Parameters
  const double minRadius = 10.0;
  const double maxRadius = 15.0;
  const double velMod = 10.0;
  // Init Coords
  double currRadius = 0.0;
  double radialVector[3] = {0.0};
  double tangVector[3] = {0.0};
  // Find the Centre Of the Domain
  double centrePoint[3] = {0.0};
  centrePoint[0] = 0.5 * (domainSizeMax[0]+domainSizeMin[0]);
  centrePoint[1] = 0.5 * (domainSizeMax[1]+domainSizeMin[1]);
  centrePoint[2] = 0.5 * (domainSizeMax[2]+domainSizeMin[2]);
  // Assign Cell Velocities
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Set to Zero
    cellPoints[loopA].concentration = 0.0;
    switch(dir){
      case kdirX:
        // Get Current Radius
        radialVector[0] = 0.0;
        radialVector[1] = cellPoints[loopA].position[1] - centrePoint[1];
        radialVector[2] = cellPoints[loopA].position[2] - centrePoint[2];
        break;
      case kdirY:
        // Get Current Radius
        radialVector[0] = cellPoints[loopA].position[0] - centrePoint[0];
        radialVector[1] = 0.0;
        radialVector[2] = cellPoints[loopA].position[2] - centrePoint[2];      
        break;
      case kdirZ:
        // Get Current Radius
        radialVector[0] = cellPoints[loopA].position[0] - centrePoint[0];
        radialVector[1] = cellPoints[loopA].position[1] - centrePoint[1];
        radialVector[2] = 0.0;
        break;
    }
    // Normalize Radial Direction
    currRadius = MRIUtils::Do3DEucNorm(radialVector);
    MRIUtils::Normalize3DVector(radialVector);
    // Eval Tangential Direction
    EvalTangentDirection(dir,radialVector,tangVector);
    // Set Velocities    
    if ((currRadius>=minRadius)&&(currRadius<=maxRadius)){
      cellPoints[loopA].velocity[0] = tangVector[0]*velMod;
      cellPoints[loopA].velocity[1] = tangVector[1]*velMod;
      cellPoints[loopA].velocity[2] = tangVector[2]*velMod;
    }
  }  
}

// ASSIGN SPHERICAL HILL VORTEX FLOW 
void MRIScan::AssignSphericalFlowSignature(MRIDirection dir){
  // Set Parameters
  const double CONST_U0 = 0.1;
  const double CONST_A = 12.0;
  // Init Local Coords
  double currentX = 0.0;
  double currentY = 0.0;
  double currentZ = 0.0;
  double localR = 0.0;
  double localZ = 0.0; 
  double axialComponentIn = 0.0;
  double radialComponentIn = 0.0;
  double axialComponentOut = 0.0;
  double radialComponentOut = 0.0;  
  // Allocate Radial and Axial Vectors
  double axialVec[3] = {0.0};
  double radialVec[3] = {0.0}; 
  // Find the Centre Of the Domain
  double centrePoint[3] = {0.0};
  centrePoint[0] = 0.5 * (domainSizeMax[0]+domainSizeMin[0]);
  centrePoint[1] = 0.5 * (domainSizeMax[1]+domainSizeMin[1]);
  centrePoint[2] = 0.5 * (domainSizeMax[2]+domainSizeMin[2]);
  // Assign Cell Velocities
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Get The Three Coordinates
    currentX = cellPoints[loopA].position[0] - centrePoint[0];
    currentY = cellPoints[loopA].position[1] - centrePoint[1];
    currentZ = cellPoints[loopA].position[2] - centrePoint[2];
    // Set a Zero Concentration
    cellPoints[loopA].concentration = 0.0;
    switch(dir){
      case kdirX:
        // Set Local Coordinates
        //localR = sqrt(currentX*currentX + currentY*currentY + currentZ*currentZ);
        localR = sqrt(currentY*currentY + currentZ*currentZ);
        localZ = currentX;
        // Build Axial Vector
        axialVec[0] = 1.0;
        axialVec[1] = 0.0;
        axialVec[2] = 0.0;
         // Build Radial Vector
        radialVec[0] = 0.0;
        radialVec[1] = currentY;
        radialVec[2] = currentZ;
        MRIUtils::Normalize3DVector(radialVec);
        break;
      case kdirY:
        // Set Local Coordinates
        localR = sqrt(currentX*currentX + currentY*currentY + currentZ*currentZ);
        localZ = currentY;                
        // Build Axial Vector
        axialVec[0] = 0.0;
        axialVec[1] = 1.0;
        axialVec[2] = 0.0;        
         // Build Radial Vector
        radialVec[0] = currentX;
        radialVec[1] = 0.0;
        radialVec[2] = currentZ;  
        MRIUtils::Normalize3DVector(radialVec);      
        break;
      case kdirZ:
        // Set Local Coordinates
        localR = sqrt(currentX*currentX + currentY*currentY + currentZ*currentZ);
        localZ = currentZ;
        // Build Axial Vector
        axialVec[0] = 0.0;
        axialVec[1] = 0.0;
        axialVec[2] = 1.0;      
         // Build Radial Vector
        radialVec[0] = currentX;
        radialVec[1] = currentY;
        radialVec[2] = 0.0;
        MRIUtils::Normalize3DVector(radialVec);
        break;
    }
    // Normalize Radial Direction
    axialComponentIn = (3.0/2.0)*CONST_U0*(1.0-((2.0*localR*localR+localZ*localZ)/(CONST_A*CONST_A)));
    radialComponentIn = (3.0/2.0)*CONST_U0*((localR*localZ)/(CONST_A*CONST_A));
    axialComponentOut = CONST_U0*(pow(((CONST_A*CONST_A)/(localR*localR+localZ*localZ)),2.5)*((2.0*localZ*localZ-localR*localR)/(2.0*CONST_A*CONST_A))-1.0);
    radialComponentOut = (3.0/2.0)*CONST_U0*((localR*localZ)/(CONST_A*CONST_A))*pow(((CONST_A*CONST_A)/(localR*localR+localZ*localZ)),2.5);
    // Set The Vector Components
    if (localR*localR+localZ*localZ<=CONST_A*CONST_A){
      cellPoints[loopA].velocity[0] = axialComponentIn * axialVec[0] + radialComponentIn * radialVec[0];
      cellPoints[loopA].velocity[1] = axialComponentIn * axialVec[1] + radialComponentIn * radialVec[1];
      cellPoints[loopA].velocity[2] = axialComponentIn * axialVec[2] + radialComponentIn * radialVec[2];
    }else{
      cellPoints[loopA].velocity[0] = 0.0;
      cellPoints[loopA].velocity[1] = 0.0;
      cellPoints[loopA].velocity[2] = 0.0;
      //cellPoints[loopA].velocity[0] = axialComponentOut * axialVec[0] + radialComponentOut * radialVec[0];
      //cellPoints[loopA].velocity[1] = axialComponentOut * axialVec[1] + radialComponentOut * radialVec[1];
      //cellPoints[loopA].velocity[2] = axialComponentOut * axialVec[2] + radialComponentOut * radialVec[2];     
    }
  }
}

// ASSIGN SPHERICAL VORTEX FLOW 
void MRIScan::AssignToroidalVortexFlowSignature(){
  // Set Parameters
  const double CONST_A = 20.0;
  const double CONST_L = 1.3*1.3;
  // Init Local Coords
  double currentX = 0.0;
  double currentY = 0.0;
  double currentZ = 0.0;
  double radius = 0.0;
  double currModulus = 0.0;
  // Allocate Radial and Axial Vectors
  double axialVec[3] = {0.0};
  double radialVec[3] = {0.0}; 
  double tangVector[3] = {0.0};
  double inclVector[3] = {0.0};
  double velVector[3] = {0.0};
  // Find the Centre Of the Domain
  double centrePoint[3] = {0.0};
  centrePoint[0] = 0.5 * (domainSizeMax[0]+domainSizeMin[0]);
  centrePoint[1] = 0.5 * (domainSizeMax[1]+domainSizeMin[1]);
  centrePoint[2] = 0.5 * (domainSizeMax[2]+domainSizeMin[2]);
  // Assign Cell Velocities
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Get The Three Coordinates
    currentX = cellPoints[loopA].position[0] - centrePoint[0];
    currentY = cellPoints[loopA].position[1] - centrePoint[1];
    currentZ = cellPoints[loopA].position[2] - centrePoint[2];
    
    // Set a Zero Concentration
    cellPoints[loopA].concentration = 0.0;
    
    // Build Axial Vector
    axialVec[0] = 1.0;
    axialVec[1] = 0.0;
    axialVec[2] = 0.0;      
    
    // Find The Projection on the YZ Plane
    radialVec[0] = 0.0;
    radialVec[1] = currentY;
    radialVec[2] = currentZ;
    MRIUtils::Normalize3DVector(radialVec);
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
        radialVec[loopB] *=  CONST_A;
    }
    
    // Find The Tangent Vector
    MRIUtils::Do3DExternalProduct(radialVec,axialVec,tangVector);
    MRIUtils::Normalize3DVector(tangVector);
    
    // Find The Inclined Vector
    inclVector[0] = currentX - radialVec[0];
    inclVector[1] = currentY - radialVec[1];
    inclVector[2] = currentZ - radialVec[2];
    radius = MRIUtils::Do3DEucNorm(inclVector);
    MRIUtils::Normalize3DVector(inclVector);
    
    // Eval Vel Vector
    MRIUtils::Do3DExternalProduct(inclVector,tangVector,velVector);
    MRIUtils::Normalize3DVector(velVector);
    
    // Eval Current Modulus
    currModulus = (8.0*radius/CONST_L)*exp(-(radius/CONST_L));
      
    // Normalize Radial Direction
    cellPoints[loopA].velocity[0] = currModulus * velVector[0];
    cellPoints[loopA].velocity[1] = currModulus * velVector[1];
    cellPoints[loopA].velocity[2] = currModulus * velVector[2];
  }
}

// ASSIGN CONSTANT FLOW
void MRIScan::AssignConstantSignature(MRIDirection dir){
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    switch(dir){
      case kdirX:
        cellPoints[loopA].velocity[0] = 1.0;
        cellPoints[loopA].velocity[1] = 0.0;
        cellPoints[loopA].velocity[2] = 0.0;
        break;
      case kdirY:
        cellPoints[loopA].velocity[0] = 0.0;
        cellPoints[loopA].velocity[1] = 1.0;
        cellPoints[loopA].velocity[2] = 0.0;
        break;
      case kdirZ:
        cellPoints[loopA].velocity[0] = 0.0;
        cellPoints[loopA].velocity[1] = 0.0;
        cellPoints[loopA].velocity[2] = 1.0;
        break;
    }
  }
}

// SET VELOCITIES TO ZERO
void MRIScan::AssignZeroVelocities(){
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    cellPoints[loopA].velocity[0] = 0.0;
    cellPoints[loopA].velocity[1] = 0.0;
    cellPoints[loopA].velocity[2] = 0.0;
  }
}

// Assign Constant Flow With Step
void MRIScan::AssignConstantFlowWithStep(){
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    if ((cellPoints[loopA].position[0]>(0.5*(domainSizeMin[0] + domainSizeMax[0])))&&
       (cellPoints[loopA].position[1]>(0.5*(domainSizeMin[1] + domainSizeMax[1])))){
      // Assign Constant Velocity
      cellPoints[loopA].concentration = 0.0;
      cellPoints[loopA].velocity[0] = 0.0;
      cellPoints[loopA].velocity[1] = 0.0;
      cellPoints[loopA].velocity[2] = 0.0;         
    }else{
      // Assign Constant Velocity
      cellPoints[loopA].concentration = 10.0;
      cellPoints[loopA].velocity[0] = 1.0;
      cellPoints[loopA].velocity[1] = 0.0;
      cellPoints[loopA].velocity[2] = 0.0;
    }
  }
}

// Assign Standard Gaussian Random Velocities on the Three Separated Components
void MRIScan::AssignRandomStandardGaussianFlow(){
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Assign Constant Velocity
    cellPoints[loopA].concentration = 10.0;
    cellPoints[loopA].velocity[0] = 1.0;
    cellPoints[loopA].velocity[1] = 0.0;
    cellPoints[loopA].velocity[2] = 0.0;
  }
}


// ASSIGN POISEILLE FLOW
void MRIScan::AssignPoiseilleSignature(MRIDirection dir){
  double currentVelocity = 0.0;
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    double currentDistance = 0.0;
    double totalDistance = 0.0;
    // Set to Zero
    cellPoints[loopA].velocity[0] = 0.0;
    cellPoints[loopA].velocity[1] = 0.0;
    cellPoints[loopA].velocity[2] = 0.0;
    cellPoints[loopA].concentration = 0.0;
    switch(dir){
      case kdirX:
        currentDistance = cellPoints[loopA].position[0];
        totalDistance = domainSizeMax[0]-domainSizeMin[0];
        break;
      case kdirY:
        currentDistance = cellPoints[loopA].position[1];
        totalDistance = domainSizeMax[1]-domainSizeMin[1];
        break;
      case kdirZ:
        currentDistance = cellPoints[loopA].position[2];
        totalDistance = domainSizeMax[2]-domainSizeMin[2];
        break;
    }
    currentVelocity = -(4.0/(totalDistance*totalDistance))*(currentDistance*currentDistance)+(4.0/totalDistance)*currentDistance;
    // Assign Velocity
    switch(dir){
      case kdirX: 
        cellPoints[loopA].velocity[1] = currentVelocity;
        break;
      case kdirY: 
        cellPoints[loopA].velocity[2] = currentVelocity;
        break;
      case kdirZ: 
        cellPoints[loopA].velocity[0] = currentVelocity;
        break;
    }
  }
}

// ASSIGN CONCENTRATIONS AND VELOCITIES
void MRIScan::AssignVelocitySignature(MRIDirection dir, MRISamples sample, double currTime){
  switch(sample){
    case kZeroVelocity:
      AssignZeroVelocities();
      break;
    case kConstantFlow:
      AssignConstantSignature(dir);
      break;    
    case kPoiseilleFlow:
      AssignPoiseilleSignature(dir);
      break;
    case kStagnationFlow:
      AssignStagnationFlowSignature(dir);
      break;
    case kCylindricalVortex:
      AssignCylindricalFlowSignature(dir);
      break;    
    case kSphericalVortex:
      AssignSphericalFlowSignature(dir);
      break;   
    case kToroidalVortex:
      AssignToroidalVortexFlowSignature();
      break;  
    case kTransientFlow:
      AssignTimeDependentPoiseilleSignature(2*3.1415,0.1,1.0e-3,currTime,1.0);
      break;
    case kConstantFlowWithStep:
      AssignConstantFlowWithStep();
      break;
  }
}

// CREATE SAMPLE FLOWS
void MRIScan::CreateSampleCase(MRISamples sampleType,
                               int sizeX, int sizeY, int sizeZ,
                               double distX, double distY, double distZ,
                               double currTime, MRIDirection dir){
  int* currentCoords = new int[kNumberOfDimensions];
  // Set Cells Totals
  cellTotals[0] = sizeX;
  cellTotals[1] = sizeY;
  cellTotals[2] = sizeZ;
  // Set Cell Lengths
  cellLength[0] = distX;
  cellLength[1] = distY;
  cellLength[2] = distZ;
  // Set Global Dimensions
  // Min
  domainSizeMin[0] = 0.0;
  domainSizeMin[1] = 0.0;
  domainSizeMin[2] = 0.0;
  // Max
  domainSizeMax[0] = (sizeX) * distX;
  domainSizeMax[1] = (sizeY) * distY;
  domainSizeMax[2] = (sizeZ) * distZ;
  // Set Total Cells
  totalCellPoints = sizeX * sizeY * sizeZ;
  // Allocate Cell Values
  cellPoints.reserve(totalCellPoints);
  // Assign Coordinates
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    MapIndexToCoords(loopA,currentCoords);
    cellPoints[loopA].position[0] = currentCoords[0] * distX;
    cellPoints[loopA].position[1] = currentCoords[1] * distY;
    cellPoints[loopA].position[2] = currentCoords[2] * distZ;
  }
  // Assign Concentrations and Velocities
  AssignVelocitySignature(dir,sampleType,currTime);
  // Find Velocity Modulus
  maxVelModule = 0.0;
  double currentMod = 0.0;
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    currentMod = sqrt((cellPoints[loopA].velocity[0])*(cellPoints[loopA].velocity[0])+
                      (cellPoints[loopA].velocity[1])*(cellPoints[loopA].velocity[1])+
                      (cellPoints[loopA].velocity[2])*(cellPoints[loopA].velocity[2]));
    if(currentMod>maxVelModule) maxVelModule = currentMod;
  }
  // Write Statistics
  std::string Msgs = WriteStatistics();
  WriteSchMessage(Msgs);
  // Deallocate
  delete [] currentCoords;
}

// ASSIGN TIME DEPENDENT FLOW
void MRIScan::AssignTimeDependentPoiseilleSignature(double omega, double radius, double viscosity, double currtime, double maxVel){
  // Eval omegaMod
  double omegaMod = ((omega*radius*radius)/viscosity);
  double relCoordX = 0.0;
  double relCoordY = 0.0;
  double relCoordZ = 0.0;
  double localRadius = 0.0;
  double normRadius = 0.0;
  double bValue = 0.0;
  // Get center point of domain
  double centrePoint[3] = {0.0};
  centrePoint[0] = 0.5 * (domainSizeMax[0] + domainSizeMin[0]);
  centrePoint[1] = 0.5 * (domainSizeMax[1] + domainSizeMin[1]);
  centrePoint[2] = 0.5 * (domainSizeMax[2] + domainSizeMin[2]);
  // Loop Through the Points
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    //Get Radius
    relCoordX = cellPoints[loopA].position[0] - centrePoint[0];
    relCoordY = cellPoints[loopA].position[1] - centrePoint[1];
    relCoordZ = cellPoints[loopA].position[2] - centrePoint[2];
    localRadius = sqrt(relCoordY*relCoordY+relCoordZ*relCoordZ);
    normRadius = (localRadius/radius);
    bValue = (1.0 - normRadius)*sqrt(omegaMod/2.0);
    // Assign Velocity
    if (normRadius<=1.0){
      if (normRadius>kMathZero){
        cellPoints[loopA].velocity[0] = maxVel*((4.0/omegaMod)*(sin(omega*currtime)-((exp(-bValue))/(sqrt(normRadius)))*sin(omega*currtime-bValue)));
      }else{
        cellPoints[loopA].velocity[0] = maxVel*((4.0/omegaMod)*(sin(omega*currtime)));
      }
      cellPoints[loopA].velocity[1] = 0.0;
      cellPoints[loopA].velocity[2] = 0.0;
    }else{
      cellPoints[loopA].velocity[0] = 0.0;
      cellPoints[loopA].velocity[1] = 0.0;
      cellPoints[loopA].velocity[2] = 0.0;      
    }
  }
}
