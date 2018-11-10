# include "mriScan.h"

using namespace std;

// STAGNATION FLOW SOLUTION
void mriScan::assignStagnationFlowSignature(mriDirection dir){
  double bConst = -1.0;
  double xCoord,yCoord,zCoord;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Set to Zero
    cells[loopA].concentration = 0.0;
    switch(dir){
      case kdirX:
        yCoord = topology->cellLocations[loopA][1] - topology->domainSizeMin[1];
        zCoord = topology->cellLocations[loopA][2] - 0.5*(topology->domainSizeMin[2]+topology->domainSizeMax[2]);
        cells[loopA].velocity[0] = 0.0;
        cells[loopA].velocity[1] = -bConst * yCoord;
        cells[loopA].velocity[2] = bConst * zCoord;
        break;
      case kdirY:
        xCoord = topology->cellLocations[loopA][0] - topology->domainSizeMin[0];
        zCoord = topology->cellLocations[loopA][2]- 0.5 * (topology->domainSizeMin[2] + topology->domainSizeMax[2]);
        cells[loopA].velocity[0] = bConst * xCoord;
        cells[loopA].velocity[1] = 0.0;
        cells[loopA].velocity[2] = -bConst * zCoord;
        break;
      case kdirZ:
        xCoord = topology->cellLocations[loopA][0] - topology->domainSizeMin[0];
        yCoord = topology->cellLocations[loopA][1] - 0.5 * (topology->domainSizeMin[1] + topology->domainSizeMax[1]);
        cells[loopA].velocity[0] = bConst * xCoord;
        cells[loopA].velocity[1] = -bConst * yCoord;
        cells[loopA].velocity[2] = 0.0;
        break;
    }
    cells[loopA].concentration = 1.0;
  }
}

// EVAL TANGENT DIRECTION
void evalTangentDirection(mriDirection dir, const mriDoubleVec& radialVector, mriDoubleVec& tangVector){
  // Create Axial Direction
  mriDoubleVec axialVec(3,0.0);
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
  mriUtils::do3DExternalProduct(axialVec,radialVector,tangVector);
  // Normalize Result
  mriUtils::normalize3DVector(tangVector);
}

// ASSIGN CYLINDRICAL VORTEX FLOW 
void mriScan::assignCylindricalFlowSignature(mriDirection dir, const mriDoubleVec& auxParams){
  // Set Parameters
  // Get Min Dimension
  double minDist = min((topology->domainSizeMax[0]-topology->domainSizeMin[0]),
                   min((topology->domainSizeMax[1]-topology->domainSizeMin[1]),
                       (topology->domainSizeMax[2]-topology->domainSizeMin[2])));
  // Set the Minimum and Maximum Radius
  const double minRadius = auxParams[0];
  const double maxRadius = auxParams[1];
  const double velMod = auxParams[2];
 // Init Coords
  double currRadius = 0.0;
  mriDoubleVec radialVector(3,0.0);
  mriDoubleVec tangVector(3,0.0);
  // Find the Centre Of the Domain
  mriDoubleVec centrePoint(3,0.0);
  centrePoint[0] = 0.5 * (topology->domainSizeMax[0] + topology->domainSizeMin[0]);
  centrePoint[1] = 0.5 * (topology->domainSizeMax[1] + topology->domainSizeMin[1]);
  centrePoint[2] = 0.5 * (topology->domainSizeMax[2] + topology->domainSizeMin[2]);
  // Assign Cell Velocities
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    switch(dir){
      case kdirX:
        // Get Current Radius
        radialVector[0] = 0.0;
        radialVector[1] = topology->cellLocations[loopA][1] - centrePoint[1];
        radialVector[2] = topology->cellLocations[loopA][2] - centrePoint[2];
        break;
      case kdirY:
        // Get Current Radius
        radialVector[0] = topology->cellLocations[loopA][0] - centrePoint[0];
        radialVector[1] = 0.0;
        radialVector[2] = topology->cellLocations[loopA][2] - centrePoint[2];      
        break;
      case kdirZ:
        // Get Current Radius
        radialVector[0] = topology->cellLocations[loopA][0] - centrePoint[0];
        radialVector[1] = topology->cellLocations[loopA][1] - centrePoint[1];
        radialVector[2] = 0.0;
        break;
    }
    // Normalize Radial Direction
    currRadius = mriUtils::do3DEucNorm(radialVector);
    mriUtils::normalize3DVector(radialVector);
    // Eval Tangential Direction
    evalTangentDirection(dir,radialVector,tangVector);
    // Set Velocities    
    if ((currRadius>=minRadius)&&(currRadius<=maxRadius)){
      cells[loopA].velocity[0] = tangVector[0]*velMod;
      cells[loopA].velocity[1] = tangVector[1]*velMod;
      cells[loopA].velocity[2] = tangVector[2]*velMod;
      cells[loopA].concentration = 1.0;
    }else{
      cells[loopA].concentration = 0.0;
    }
  }  
}

// ASSIGN SPHERICAL HILL VORTEX FLOW 
void mriScan::assignSphericalFlowSignature(mriDirection dir, const mriDoubleVec& auxParams){
  // Set Parameters
  double minDist = min(topology->domainSizeMax[0]-topology->domainSizeMin[0],
                   min(topology->domainSizeMax[1]-topology->domainSizeMin[1],
                       topology->domainSizeMax[2]-topology->domainSizeMin[2]));
  const double CONST_U0 = auxParams[0];
  const double CONST_A = auxParams[1];
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
  mriDoubleVec axialVec(3,0.0);
  mriDoubleVec radialVec(3,0.0);
  // Find the Centre Of the Domain
  mriDoubleVec centrePoint(3,0.0);
  centrePoint[0] = 0.5 * (topology->domainSizeMax[0]+topology->domainSizeMin[0]);
  centrePoint[1] = 0.5 * (topology->domainSizeMax[1]+topology->domainSizeMin[1]);
  centrePoint[2] = 0.5 * (topology->domainSizeMax[2]+topology->domainSizeMin[2]);
  // Assign Cell Velocities
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Get The Three Coordinates
    currentX = topology->cellLocations[loopA][0] - centrePoint[0];
    currentY = topology->cellLocations[loopA][1] - centrePoint[1];
    currentZ = topology->cellLocations[loopA][2] - centrePoint[2];
    switch(dir){
      case kdirX:
        // Set Local Coordinates
        //localR = sqrt(currentX*currentX + currentY*currentY + currentZ*currentZ);
        localR = sqrt(currentY * currentY + currentZ * currentZ);
        localZ = currentX;
        // Build Axial Vector
        axialVec[0] = 1.0;
        axialVec[1] = 0.0;
        axialVec[2] = 0.0;
         // Build Radial Vector
        radialVec[0] = 0.0;
        radialVec[1] = currentY;
        radialVec[2] = currentZ;
        mriUtils::normalize3DVector(radialVec);
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
        mriUtils::normalize3DVector(radialVec);      
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
        mriUtils::normalize3DVector(radialVec);
        break;
    }
    // Assign Concentration
    cells[loopA].concentration = 1.0;
    // Normalize Radial Direction
    axialComponentIn = (3.0/2.0)*CONST_U0*(1.0-((2.0*localR*localR+localZ*localZ)/(CONST_A*CONST_A)));
    radialComponentIn = (3.0/2.0)*CONST_U0*((localR*localZ)/(CONST_A*CONST_A));
    axialComponentOut = CONST_U0*(pow(((CONST_A*CONST_A)/(localR*localR+localZ*localZ)),2.5)*((2.0*localZ*localZ-localR*localR)/(2.0*CONST_A*CONST_A))-1.0);
    radialComponentOut = (3.0/2.0)*CONST_U0*((localR*localZ)/(CONST_A*CONST_A))*pow(((CONST_A*CONST_A)/(localR*localR+localZ*localZ)),2.5);
    // Set The Vector Components
    if (localR*localR+localZ*localZ<=CONST_A*CONST_A){
      // Set a Zero Concentration
      cells[loopA].concentration = 1.0;
      // Set Velocities      
      cells[loopA].velocity[0] = axialComponentIn * axialVec[0] + radialComponentIn * radialVec[0];
      cells[loopA].velocity[1] = axialComponentIn * axialVec[1] + radialComponentIn * radialVec[1];
      cells[loopA].velocity[2] = axialComponentIn * axialVec[2] + radialComponentIn * radialVec[2];
    }else{
      // Set a Zero Concentration
      cells[loopA].concentration = 0.0;
      // Set Velocities
      cells[loopA].velocity[0] = 0.0;
      cells[loopA].velocity[1] = 0.0;
      cells[loopA].velocity[2] = 0.0;
      //cellPoints[loopA].velocity[0] = axialComponentOut * axialVec[0] + radialComponentOut * radialVec[0];
      //cellPoints[loopA].velocity[1] = axialComponentOut * axialVec[1] + radialComponentOut * radialVec[1];
      //cellPoints[loopA].velocity[2] = axialComponentOut * axialVec[2] + radialComponentOut * radialVec[2];     
    }
  }
}

// ASSIGN SPHERICAL VORTEX FLOW 
void mriScan::assignToroidalVortexFlowSignature(const mriDoubleVec& auxParams){
  // Set Parameters
  const double CONST_A = auxParams[0];
  const double CONST_L = auxParams[1];
  // Init Local Coords
  double currentX = 0.0;
  double currentY = 0.0;
  double currentZ = 0.0;
  double radius = 0.0;
  double currModulus = 0.0;
  // Allocate Radial and Axial Vectors
  mriDoubleVec axialVec(3,0.0);
  mriDoubleVec radialVec(3,0.0); 
  mriDoubleVec tangVector(3,0.0);
  mriDoubleVec inclVector(3,0.0);
  mriDoubleVec velVector(3,0.0);
  // Find the Centre Of the Domain
  mriDoubleVec centrePoint(3,0.0);
  centrePoint[0] = 0.5 * (topology->domainSizeMax[0]+topology->domainSizeMin[0]);
  centrePoint[1] = 0.5 * (topology->domainSizeMax[1]+topology->domainSizeMin[1]);
  centrePoint[2] = 0.5 * (topology->domainSizeMax[2]+topology->domainSizeMin[2]);
  // Assign Cell Velocities
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Get The Three Coordinates
    currentX = topology->cellLocations[loopA][0] - centrePoint[0];
    currentY = topology->cellLocations[loopA][1] - centrePoint[1];
    currentZ = topology->cellLocations[loopA][2] - centrePoint[2];
    
    // Set a Zero Concentration
    cells[loopA].concentration = 1.0;
    
    // Build Axial Vector
    axialVec[0] = 1.0;
    axialVec[1] = 0.0;
    axialVec[2] = 0.0;      
    
    // Find The Projection on the YZ Plane
    radialVec[0] = 0.0;
    radialVec[1] = currentY;
    radialVec[2] = currentZ;
    mriUtils::normalize3DVector(radialVec);
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
        radialVec[loopB] *=  CONST_A;
    }
    
    // Find The Tangent Vector
    mriUtils::do3DExternalProduct(radialVec,axialVec,tangVector);
    mriUtils::normalize3DVector(tangVector);
    
    // Find The Inclined Vector
    inclVector[0] = currentX - radialVec[0];
    inclVector[1] = currentY - radialVec[1];
    inclVector[2] = currentZ - radialVec[2];
    radius = mriUtils::do3DEucNorm(inclVector);
    mriUtils::normalize3DVector(inclVector);
    
    // Eval Vel Vector
    mriUtils::do3DExternalProduct(inclVector,tangVector,velVector);
    mriUtils::normalize3DVector(velVector);
    
    // Eval Current Modulus
    currModulus = (8.0*radius/CONST_L)*exp(-(radius/CONST_L));
      
    // Normalize Radial Direction
    cells[loopA].velocity[0] = currModulus * velVector[0];
    cells[loopA].velocity[1] = currModulus * velVector[1];
    cells[loopA].velocity[2] = currModulus * velVector[2];
  }
}

// ASSIGN CONSTANT FLOW
void mriScan::assignConstantSignature(mriDirection dir){
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    switch(dir){
      case kdirX:
        cells[loopA].velocity[0] = 1.0;
        cells[loopA].velocity[1] = 0.0;
        cells[loopA].velocity[2] = 0.0;
        break;
      case kdirY:
        cells[loopA].velocity[0] = 0.0;
        cells[loopA].velocity[1] = 1.0;
        cells[loopA].velocity[2] = 0.0;
        break;
      case kdirZ:
        cells[loopA].velocity[0] = 0.0;
        cells[loopA].velocity[1] = 0.0;
        cells[loopA].velocity[2] = 1.0;
        break;
    }
  }
}

// ====================
// TAYLOR FLOW VORTEX
// ====================
void mriScan::assignTaylorVortexSignature(mriDirection dir){
  
  mriDoubleVec centrePoint(3,0.0);

  centrePoint[0] = 0.5 * (topology->domainSizeMax[0] + topology->domainSizeMin[0]);
  centrePoint[1] = 0.5 * (topology->domainSizeMax[1] + topology->domainSizeMin[1]);
  centrePoint[2] = 0.5 * (topology->domainSizeMax[2] + topology->domainSizeMin[2]);
  
  // Loop on cell
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    switch(dir){
      case kdirX:
        cells[loopA].velocity[0] = 0.0;
        cells[loopA].velocity[1] = topology->cellLocations[loopA][2]-centrePoint[2];
        cells[loopA].velocity[2] = -(topology->cellLocations[loopA][1]-centrePoint[1]);
        break;
      case kdirY:
        cells[loopA].velocity[0] = topology->cellLocations[loopA][2]-centrePoint[2];
        cells[loopA].velocity[1] = 0.0;
        cells[loopA].velocity[2] = -(topology->cellLocations[loopA][0]-centrePoint[0]);
        break;
      case kdirZ:
        cells[loopA].velocity[0] = topology->cellLocations[loopA][1]-centrePoint[1];
        cells[loopA].velocity[1] = -(topology->cellLocations[loopA][0]-centrePoint[0]);
        cells[loopA].velocity[2] = 0.0;
        break;
    }
    cells[loopA].concentration = 1.0;
  }
}


// SET VELOCITIES TO ZERO
void mriScan::assignZeroVelocities(){
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    cells[loopA].velocity[0] = 0.0;
    cells[loopA].velocity[1] = 0.0;
    cells[loopA].velocity[2] = 0.0;
  }
}

// Assign Constant Flow With Step
void mriScan::assignConstantFlowWithStep(){
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    if ((topology->cellLocations[loopA][0] > (0.5*(topology->domainSizeMin[0] + topology->domainSizeMax[0])))&&
       (topology->cellLocations[loopA][1] > (0.5*(topology->domainSizeMin[1] + topology->domainSizeMax[1])))){
      // Assign Constant Velocity
      cells[loopA].concentration = 0.0;
      cells[loopA].velocity[0] = 0.0;
      cells[loopA].velocity[1] = 0.0;
      cells[loopA].velocity[2] = 0.0;         
    }else{
      // Assign Constant Velocity
      cells[loopA].concentration = 10.0;
      cells[loopA].velocity[0] = 1.0;
      cells[loopA].velocity[1] = 0.0;
      cells[loopA].velocity[2] = 0.0;
    }
  }
}

// Assign Standard Gaussian Random Velocities on the Three Separated Components
void mriScan::assignRandomStandardGaussianFlow(){
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Assign Constant Velocity
    cells[loopA].concentration = 10.0;
    cells[loopA].velocity[0] = 1.0;
    cells[loopA].velocity[1] = 0.0;
    cells[loopA].velocity[2] = 0.0;
  }
}

// =====================
// ASSIGN POISEILLE FLOW
// ===================== 
void mriScan::assignPoiseuilleSignature(mriDirection dir, const mriDoubleVec& auxParams){
  double totalDistance   = auxParams[0];
  double peakVel         = auxParams[1];
  double currentVelocity = 0.0;
  double conc            = 0.0;
  // SET CENTER POINT
  double centerPoint[3];
  centerPoint[0] = 0.5*(topology->domainSizeMax[0] + topology->domainSizeMin[0]);
  centerPoint[1] = 0.5*(topology->domainSizeMax[1] + topology->domainSizeMin[1]);
  centerPoint[2] = 0.5*(topology->domainSizeMax[2] + topology->domainSizeMin[2]);
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    double currentDistance = 0.0;
    // Set to Zero
    cells[loopA].velocity[0]   = 0.0;
    cells[loopA].velocity[1]   = 0.0;
    cells[loopA].velocity[2]   = 0.0;
    cells[loopA].concentration = 0.0;
    switch(dir){
      case kdirX:
        currentDistance = sqrt((topology->cellLocations[loopA][1] - centerPoint[1])*(topology->cellLocations[loopA][1] - centerPoint[1]) +
                               (topology->cellLocations[loopA][2] - centerPoint[2])*(topology->cellLocations[loopA][2] - centerPoint[2]));
        break;
      case kdirY:
        currentDistance = sqrt((topology->cellLocations[loopA][0] - centerPoint[0])*(topology->cellLocations[loopA][0] - centerPoint[0]) +
                               (topology->cellLocations[loopA][2] - centerPoint[2])*(topology->cellLocations[loopA][2] - centerPoint[2]));
        break;
      case kdirZ:
        currentDistance = sqrt((topology->cellLocations[loopA][0] - centerPoint[0])*(topology->cellLocations[loopA][0] - centerPoint[0]) +
                               (topology->cellLocations[loopA][1] - centerPoint[1])*(topology->cellLocations[loopA][1] - centerPoint[1]));
        break;
    }
    // Apply a threshold
    if(currentDistance<totalDistance){
      currentVelocity = -(peakVel/(totalDistance*totalDistance))*(currentDistance*currentDistance) + peakVel;
      conc = 1.0;
    }else{
      currentVelocity = 0.0;
      conc = 1.0e-10;
    }
    // Assign Velocity
    cells[loopA].concentration = conc;
    switch(dir){
      case kdirX: 
        cells[loopA].velocity[0] = currentVelocity;
        break;
      case kdirY: 
        cells[loopA].velocity[1] = currentVelocity;
        break;
      case kdirZ: 
        cells[loopA].velocity[2] = currentVelocity;
        break;
    }
  }
}

// ASSIGN CONCENTRATIONS AND VELOCITIES
void mriScan::assignVelocitySignature(mriDirection dir, mriTemplateType sample, double currTime, const mriDoubleVec& auxParams){
  switch(sample){
    case kZeroVelocity:
      assignZeroVelocities();
      break;
    case kConstantFlow:
      assignConstantSignature(dir);
      break;    
    case kPoiseuilleFlow:
      assignPoiseuilleSignature(dir,auxParams);
      break;
    case kStagnationFlow:
      assignStagnationFlowSignature(dir);
      break;
    case kCylindricalVortex:
      assignCylindricalFlowSignature(dir,auxParams);
      break;    
    case kSphericalVortex:
      assignSphericalFlowSignature(dir,auxParams);
      break;   
    case kToroidalVortex:
      assignToroidalVortexFlowSignature(auxParams);
      break;  
    case kTransientFlow:
      assignTimeDependentPoiseilleSignature(2*3.1415,8.0,1.0e-3,currTime,1.0);
      break;
    case kConstantFlowWithStep:
      assignConstantFlowWithStep();
      break;
    case kTaylorVortex:
      assignTaylorVortexSignature(dir);
      break;
  }
}

// =====================
// CREATE TEMPLATE FLOWS
// =====================
void mriScan::createFromTemplate(mriTemplateType sampleType,const mriDoubleVec& params){

  // Store Parameter Values
  int sizeX       = int(params[0]);
  int sizeY       = int(params[1]);
  int sizeZ       = int(params[2]);
  double distX    = params[3];
  double distY    = params[4];
  double distZ    = params[5];
  double currTime = params[6];  

  // Slice the additional parameters if any
  mriDoubleVec::const_iterator first;
  mriDoubleVec::const_iterator last;
  if(params.size()>8){
    first = params.begin() + 8;
    last  = params.end();
  }else{
    first = params.end();
    last  = params.end();    
  }
  mriDoubleVec auxParams(first,last);

  // Template Orientation
  int dir = 0;
  int direction = int(params[7]);
  if(direction == 0){
    dir = kdirX;
  }else if(direction == 1){
    dir = kdirY;
  }else if(direction == 2){
    dir = kdirZ;
  }else{
    throw mriException("ERROR: Invalid template direction in CreateSampleCase.\n");
  }

  //Allocate Based on the Topology
  cells.clear();

  mriCell cell;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    cells.push_back(cell);
  }

  // Assign Concentrations and Velocities
  assignVelocitySignature(dir,sampleType,currTime,auxParams);

  // Find Velocity Modulus
  maxVelModule = 0.0;
  double currentMod = 0.0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    currentMod = sqrt((cells[loopA].velocity[0])*(cells[loopA].velocity[0])+
                      (cells[loopA].velocity[1])*(cells[loopA].velocity[1])+
                      (cells[loopA].velocity[2])*(cells[loopA].velocity[2]));
    if(currentMod>maxVelModule) maxVelModule = currentMod;
  }
}

// ASSIGN TIME DEPENDENT FLOW
void mriScan::assignTimeDependentPoiseilleSignature(double omega, double radius, double viscosity, double currtime, double maxVel){
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
  centrePoint[0] = 0.5 * (topology->domainSizeMax[0] + topology->domainSizeMin[0]);
  centrePoint[1] = 0.5 * (topology->domainSizeMax[1] + topology->domainSizeMin[1]);
  centrePoint[2] = 0.5 * (topology->domainSizeMax[2] + topology->domainSizeMin[2]);
  // Loop Through the Points
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    //Get Radius
    relCoordX = topology->cellLocations[loopA][0] - centrePoint[0];
    relCoordY = topology->cellLocations[loopA][1] - centrePoint[1];
    relCoordZ = topology->cellLocations[loopA][2] - centrePoint[2];
    localRadius = sqrt(relCoordY*relCoordY+relCoordZ*relCoordZ);
    normRadius = (localRadius/radius);
    bValue = (1.0 - normRadius)*sqrt(omegaMod/2.0);
    // Assign Velocity
    double peakVel = 4.0;
    if (normRadius<=1.0){
      if (normRadius>kMathZero){
        cells[loopA].velocity[0] = maxVel*((peakVel/omegaMod)*(sin(omega*currtime)-((exp(-bValue))/(sqrt(normRadius)))*sin(omega*currtime-bValue)));
      }else{
        cells[loopA].velocity[0] = maxVel*((peakVel/omegaMod)*(sin(omega*currtime)));
      }
      cells[loopA].velocity[1] = 0.0;
      cells[loopA].velocity[2] = 0.0;
      cells[loopA].concentration = 1.0;
    }else{
      cells[loopA].velocity[0] = 0.0;
      cells[loopA].velocity[1] = 0.0;
      cells[loopA].velocity[2] = 0.0;
      cells[loopA].concentration = 0.0;
    }
  }
}
