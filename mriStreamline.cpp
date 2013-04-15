#include <math.h>

#include "mriScan.h"
#include "mriStreamlineOptions.h"
#include "mriStreamline.h"
#include "mriUtils.h"

MRIStreamline::MRIStreamline()
{
}

MRIStreamline::~MRIStreamline()
{
}

// FIND SL INTERSECTION WITH FIXED X,Y,Z COORD
bool MRIStreamline::EvalSLIntersection(MRIDirection dir, double fixedCoord, double* &intCoords){
  bool found = false;
  int count = 0;
  bool isBiggerThenMin = false;  
  bool isSmallerThanMax = false;
  double currentCoord = 0.0;
  double nextCoord = 0.0;
  while ((!found)&&(count<(totalPoints-1))){
    switch(dir){
      case kdirX:
        currentCoord = xCoords[count];
        nextCoord = xCoords[count+1];
        break;
      case kdirY:
        currentCoord = yCoords[count];
        nextCoord = yCoords[count+1];
        break;
      case kdirZ:
        currentCoord = zCoords[count];
        nextCoord = zCoords[count+1];
        break;      
    } 
    // Check If limits are met
    if (fabs(currentCoord)<kMathZero){
      isBiggerThenMin = (fixedCoord>=currentCoord);
    }else{
      isBiggerThenMin = (fixedCoord>currentCoord);
    }
    isSmallerThanMax = (fixedCoord<=nextCoord);
    
    // Check If Found
    found = ((isBiggerThenMin)&&(isSmallerThanMax));
    
    // Update 
    count++;
  }
  // Exit If Not Itersection has Been Found
  if(!found){
    return found;
  }
  // Find Ratio
  double ratio = (fixedCoord-currentCoord)/(nextCoord-currentCoord);
  // Find Versor
  double currentVersor[3] = {0.0};
  currentVersor[0] = xCoords[count+1]-xCoords[count];
  currentVersor[1] = yCoords[count+1]-yCoords[count];
  currentVersor[2] = zCoords[count+1]-zCoords[count];
  MRIUtils::Normalize3DVector(currentVersor);
  // Write Result
  intCoords[0] = xCoords[count] + ratio * currentVersor[0];
  intCoords[1] = yCoords[count] + ratio * currentVersor[1];
  return found;
}

// Perform Three-linear Interpolation
void MRIScan::PerformVelocityLinearInterpolation(double* coords, int currentCell, 
                                                 int xCell, int yCell, int zCell,
                                                 double* &velocity){
                                          
  double pointRelX = coords[0]-cellPoints[currentCell].position[0];
  double pointRelY = coords[1]-cellPoints[currentCell].position[1];
  double pointRelZ = coords[2]-cellPoints[currentCell].position[2];
  double distX = cellPoints[xCell].position[0] - cellPoints[currentCell].position[0];
  double distY = cellPoints[yCell].position[1] - cellPoints[currentCell].position[1];
  double distZ = cellPoints[zCell].position[2] - cellPoints[currentCell].position[2];
  // Init Ratios
  double ratioX = 0.0;
  double ratioY = 0.0;
  double ratioZ = 0.0;  
  // X Cell
  if (xCell == currentCell){
    velocity[0] = cellPoints[currentCell].velocity[0];
  }else{
    ratioX = (pointRelX/distX);
    velocity[0] = (1.0-ratioX)*cellPoints[currentCell].velocity[0]+(ratioX)*cellPoints[xCell].velocity[0];
  }
  // Y Cell
  if (yCell == currentCell){
    velocity[1] = cellPoints[currentCell].velocity[1];
  }else{
    ratioY = (pointRelY/distY);
    velocity[1] = (1.0-ratioY)*cellPoints[currentCell].velocity[1]+(ratioY)*cellPoints[yCell].velocity[1];
  }
  // Z Cell
  if (zCell == currentCell){
    velocity[2] = cellPoints[currentCell].velocity[2];
  }else{
    ratioZ = (pointRelZ/distZ);
    velocity[2] = (1.0-ratioZ)*cellPoints[currentCell].velocity[2]+(ratioZ)*cellPoints[zCell].velocity[2];
  }
}

// Get Point Velocity
void MRIScan::GetPointVelocity(double xCoord, double yCoord, double zCoord,
                               double* &pointVel){

  // Get Coords
  double coords[3] = {0.0};
  int* currentCellIndexes = new int[3];
  coords[0] = xCoord;
  coords[1] = yCoord;
  coords[2] = zCoord;
  int xCell = 0;
  int yCell = 0;
  int zCell = 0;
  // Cell Indexes
  double xCellIndex[3] = {0.0};
  double yCellIndex[3] = {0.0};
  double zCellIndex[3] = {0.0};  
  // Find The Cells It Belongs To
  int currentCell = GetCellNumber(coords);
  // Find the Index of the Current Cell
  MapIndexToCoords(currentCell,currentCellIndexes);
  // Eval the Differences
  double xDiff = coords[0] - cellPoints[currentCell].position[0];
  double yDiff = coords[1] - cellPoints[currentCell].position[1];
  double zDiff = coords[2] - cellPoints[currentCell].position[2];
  // Normalize Differences to +1/-1
  if (fabs(xDiff)>kMathZero) xDiff = (xDiff/fabs(xDiff));
  else xDiff = 0.0;
  if (fabs(yDiff)>kMathZero) yDiff = (yDiff/fabs(yDiff));
  else yDiff = 0.0;
  if (fabs(zDiff)>kMathZero) zDiff = (zDiff/fabs(zDiff));
  else zDiff = 0.0;
  // Find Neighbour indexes
  // X
  xCellIndex[0] = currentCellIndexes[0]+(int)(xDiff);
  xCellIndex[1] = currentCellIndexes[1];
  xCellIndex[2] = currentCellIndexes[2];
  if ((xCellIndex[0]<0)||(xCellIndex[0]>cellTotals[0]-1)) xCell = currentCell;
  else xCell = MapCoordsToIndex(xCellIndex[0],xCellIndex[1],xCellIndex[2]);
  // Y
  yCellIndex[0] = currentCellIndexes[0];
  yCellIndex[1] = currentCellIndexes[1]+ (int)(yDiff);
  yCellIndex[2] = currentCellIndexes[2];
  if ((yCellIndex[1]<0)||(yCellIndex[1]>cellTotals[1]-1)) yCell = currentCell;
  else yCell = MapCoordsToIndex(yCellIndex[0],yCellIndex[1],yCellIndex[2]);
  // Z
  zCellIndex[0] = currentCellIndexes[0];
  zCellIndex[1] = currentCellIndexes[1];
  zCellIndex[2] = currentCellIndexes[2]+(int)(zDiff);
  if ((zCellIndex[2]<0)||(zCellIndex[2]>cellTotals[2]-1)) zCell = currentCell;
  else zCell = MapCoordsToIndex(zCellIndex[0],zCellIndex[1],zCellIndex[2]);
  // Get Velocity
  PerformVelocityLinearInterpolation(coords,currentCell,xCell,yCell,zCell,pointVel);
}

// Eval Single Streamline
void MRIScan::EvalSingleStreamLine(double* start, MRIStreamlineOptions &options, MRIStreamline* &sL){
  const double lengthTol = 1.0e-3;
  double newPoint[3] = {0.0};
  double* eulerVel = new double[3];
  double* otherVel = new double[3];
  // Eval Total
  int numPoints = (MRIUtils::round(options.totalT/options.deltaT)+1);
  // Allocate
  double* localxCoords = new double[numPoints];
  double* localyCoords = new double[numPoints];
  double* localzCoords = new double[numPoints];
  // Characteristic Length of the Domain
  double referenceLength = maxVelModule*options.deltaT;
  // Set The First Point To The Start Position
  localxCoords[0] = start[0];
  localyCoords[0] = start[1];
  localzCoords[0] = start[2];
  // Loop On all Time Steps
  bool finished = false;
  int count = 0;
  while ((!finished)&&(count<numPoints-1)){
    // USE RUNGE KUTTA SECOND ORDER SCHEME

    // Get the Velocity At the Current Location
    GetPointVelocity(localxCoords[count],localyCoords[count],localzCoords[count],eulerVel);

    // Eval New Location Based on Extracted Velocities
    newPoint[0] = localxCoords[count]+options.deltaT*eulerVel[0];
    newPoint[1] = localyCoords[count]+options.deltaT*eulerVel[1];
    newPoint[2] = localzCoords[count]+options.deltaT*eulerVel[2];

    // Check If outside the Domain
    if ((newPoint[0]<domainSizeMin[0])||(newPoint[0]>domainSizeMax[0])) finished = true;
    if ((newPoint[1]<domainSizeMin[1])||(newPoint[1]>domainSizeMax[1])) finished = true;
    if ((newPoint[2]<domainSizeMin[2])||(newPoint[2]>domainSizeMax[2])) finished = true;

    // Eval Velocities at New Location
    GetPointVelocity(newPoint[0],newPoint[1],newPoint[2],otherVel);

    // Eval New Point
    localxCoords[count+1] = localxCoords[count]+options.deltaT*(0.5*(eulerVel[0]+otherVel[0]));
    localyCoords[count+1] = localyCoords[count]+options.deltaT*(0.5*(eulerVel[1]+otherVel[1]));
    localzCoords[count+1] = localzCoords[count]+options.deltaT*(0.5*(eulerVel[2]+otherVel[2]));

    // Check If outside the Domain
    if ((localxCoords[count]<domainSizeMin[0])||(localxCoords[count]>domainSizeMax[0])) finished = true;
    if ((localyCoords[count]<domainSizeMin[1])||(localyCoords[count]>domainSizeMax[1])) finished = true;
    if ((localzCoords[count]<domainSizeMin[2])||(localzCoords[count]>domainSizeMax[2])) finished = true;

    // Check If Increment Is Zero
    if ((((localxCoords[count+1]-localxCoords[count])/referenceLength)<lengthTol)&&
        (((localyCoords[count+1]-localyCoords[count])/referenceLength)<lengthTol)&&
        (((localzCoords[count+1]-localzCoords[count])/referenceLength)<lengthTol)) finished = true;

    // Update Counter
    count++;
  }
  // Add positions to Streamlines
  sL->totalPoints = count;
  for(int loopA=0;loopA<count;loopA++){
    // Add x
    sL->xCoords.push_back(localxCoords[loopA]);
    // Add y
    sL->yCoords.push_back(localyCoords[loopA]);
    // Add z
    sL->zCoords.push_back(localzCoords[loopA]);
  }
  // Free Vectors
  delete [] localxCoords;
  delete [] localyCoords;
  delete [] localzCoords;  
}

// Recover MaxMin Coords
void MRIScan::RecoverPlaneMaxMinCoords(MRIPlane plane, double* minCoords, double* maxCoords){
  switch(plane){
    case kPlaneXY:
      minCoords[0] = domainSizeMin[0];
      minCoords[1] = domainSizeMin[1];
      minCoords[2] = domainSizeMin[2];
      maxCoords[0] = domainSizeMax[0];
      maxCoords[1] = domainSizeMax[1];
      maxCoords[2] = domainSizeMax[2];
      break;
    case kPlaneYZ:
      minCoords[0] = domainSizeMin[1];
      minCoords[1] = domainSizeMin[2];
      minCoords[2] = domainSizeMin[0];
      maxCoords[0] = domainSizeMax[1];
      maxCoords[1] = domainSizeMax[2];
      maxCoords[2] = domainSizeMax[0];
      break;
    case kPlaneZX:
      minCoords[0] = domainSizeMin[2];
      minCoords[1] = domainSizeMin[0];
      minCoords[2] = domainSizeMin[1];
      maxCoords[0] = domainSizeMax[2];
      maxCoords[1] = domainSizeMax[0];
      maxCoords[2] = domainSizeMax[1];
      break;
  }
}

// COMPUTE STREAMLINES
void MRIScan::ComputeStreamlines(MRIStreamlineOptions &options, std::vector<MRIStreamline*> &streamlines){
  // Init
  double minCoords[3] = {0.0};
  double maxCoords[3] = {0.0};
  double minWindow[2] = {0.0};
  double maxWindow[2] = {0.0};
  double start[3] = {0.0};
  double firstCoord = 0.0;
  double secondCoord = 0.0;
  double thirdCoord = 0.0;
  // Compute Streamlines
  //streamlines.reserve(options.gridTotals[0]*options.gridTotals[1]+1);
  // Get Min and Max Coords
  RecoverPlaneMaxMinCoords(options.planeSlice,minCoords,maxCoords);
  // Get Max and Min Window
  minWindow[0] = minCoords[0] + options.minCoordFactor[0] * (maxCoords[0]-minCoords[0]);
  minWindow[1] = minCoords[1] + options.minCoordFactor[1] * (maxCoords[1]-minCoords[1]);
  maxWindow[0] = maxCoords[0] - options.minCoordFactor[0] * (maxCoords[0]-minCoords[0]);
  maxWindow[1] = maxCoords[1] - options.minCoordFactor[1] * (maxCoords[1]-minCoords[1]);
  double distance3 = (maxCoords[2] - minCoords[2]);
  // Initialize Streamline
  int currentSL = 0;
  for(int loopA=0;loopA<options.gridTotals[0];loopA++){
    for(int loopB=0;loopB<options.gridTotals[1];loopB++){
      // Increment SL Number
      currentSL++;
      // Eval The Starting Point of the Current Streamline
      firstCoord =  minWindow[0] +(maxWindow[0]-minWindow[0])*((loopA-1)/double(options.gridTotals[0]-1));
      secondCoord = minWindow[1] +(maxWindow[1]-minWindow[1])*((loopB-1)/double(options.gridTotals[1]-1));
      thirdCoord =  minCoords[2] + options.distanceFactor * distance3;
      // Set the Right Plane
      switch(options.planeSlice){
        case kPlaneXY:
          start[0] = firstCoord;
          start[1] = secondCoord;
          start[2] = thirdCoord;
          break;
        case kPlaneYZ:
          start[0] = thirdCoord;
          start[1] = firstCoord;
          start[2] = secondCoord;
          break;
        case kPlaneZX:
          start[0] = secondCoord;
          start[1] = thirdCoord;
          start[2] = firstCoord;
          break;
      }
      // Eval Streamline  
      MRIStreamline* newSL = new MRIStreamline;
      EvalSingleStreamLine(start,options,newSL);
      streamlines.push_back(newSL);
    }
  }
}

// Eval Arrival Point Statistics
void MRIScan::EvalSLArrivalPointDistribution(int totalSL, std::vector<MRIStreamline*> &streamlines,
                                             MRIDirection dir, double minCoord, double maxCoord, 
                                             int totalSlices, std::vector<double> &sliceCenter, std::vector<double> &sliceNormArrivals){
  // Init
  double sliceMin = 0.0;
  double sliceMax = 0.0;
  bool biggerThanMin = false;
  bool smallerThanMax = false;
  double lastCoord = 0.0;
  // Set Margin
  double refLength = maxCoord - minCoord;
  double margin = 0.1 * refLength;
  // Allocate
  sliceCenter.reserve(totalSlices);
  // Initialize Slice Centers
  sliceCenter[0] = minCoord + (refLength/double(2.0*totalSlices));
  for(int loopA=1;loopA<totalSlices;loopA++){
    sliceCenter[loopA] = sliceCenter[loopA-1]+(refLength/double(totalSlices));
  }
  // Init Slice Norm Arrivals
  sliceNormArrivals.reserve(totalSlices);
  for(int loopA=0;loopA<totalSlices;loopA++){ 
    sliceNormArrivals[loopA] = 0.0;
  }
  for(int loopA=0;loopA<totalSL;loopA++){
    // Get Last Coord
    switch(dir){
      case kdirX:
        lastCoord = streamlines[loopA]->xCoords[streamlines[loopA]->totalPoints-1];
        break;
      case kdirY:
        lastCoord = streamlines[loopA]->yCoords[streamlines[loopA]->totalPoints-1];
        break;
      case kdirZ:
        lastCoord = streamlines[loopA]->zCoords[streamlines[loopA]->totalPoints-1];
        break;      
    }
    for(int loopB=0;loopB<totalSlices;loopB++){
      // Count Min and Max Slices
      sliceMin = sliceCenter[loopB] - (refLength/double(2.0*totalSlices));
      sliceMax = sliceCenter[loopB] + (refLength/double(2.0*totalSlices));
      // Check Value Within Bounds
      // Min Value
      if (fabs(sliceMin-minCoord)>kMathZero) {
        biggerThanMin = (lastCoord>sliceMin);
      }else{
        biggerThanMin = (lastCoord>=sliceMin-margin);
      }
      // Max Value
      if (fabs(sliceMax-maxCoord)>kMathZero){
        smallerThanMax = (lastCoord<sliceMax);
      }else{
        smallerThanMax = (lastCoord<=sliceMax+margin);
      }
      //  Add to Counter
      if ((biggerThanMin)&&(smallerThanMax)){
        sliceNormArrivals[loopB] = sliceNormArrivals[loopB] + 1.0;
      }
    }
  }
  // Init Counter
  int totalCollected = 0;
  for(int loopA=0;loopA<totalSlices;loopA++){
    totalCollected = totalCollected + (int)(sliceNormArrivals[loopA]);
  }
  // Check Number of Collected Streamlines
  if (totalCollected != totalSL){
    printf("Error: Some Streamlines are Missing!!!\n");
  }
  // Make Cumulative
  for(int loopA=0;loopA<totalSlices;loopA++){
    sliceNormArrivals[loopA] = sliceNormArrivals[loopA] + sliceNormArrivals[loopA-1];
  }
}

// Eval Transverse Diffusion Statistics
void EvalSLTransverseDiffusionWithTime(int totalSL, std::vector<MRIStreamline> streamlines, 
                                       MRIDirection dir, double minCoord, double maxCoord,
                                       MRIStreamlineOptions options,
                                       std::vector<double> &time, std::vector<double> &crossDeviations){
  double currentAv = 0.0;
  double currentSqrAv = 0.0;
  int count = 0;
  double lastCoord = 0.0;
  double xDiff = 0.0;
  double yDiff = 0.0;
  double dist = 0.0;
  double value = 0.0;
  // Init
  double refLength = maxCoord - minCoord;
  int totalTimeSteps = (MRIUtils::round(options.totalT/double(options.deltaT))+1);
  // Allocate and Initialize
  time.reserve(totalTimeSteps);
  crossDeviations.reserve(totalTimeSteps);
  time[0] = 0.0;
  for(int loopA=1;loopA<totalTimeSteps;loopA++){
    time[loopA] = time[loopA-1] + options.deltaT;
  }
  // Loop On Time Steps
  for(int loopA=0;loopA<totalTimeSteps;loopA++){
    currentAv = 0.0;
    currentSqrAv = 0.0;
    count = 0;
    for(int loopB=0;loopB<totalSL;loopB++){
      switch(dir){
        case kdirX: 
          lastCoord = streamlines[loopB].xCoords[streamlines[loopB].totalPoints-1];
          break;
        case kdirY: 
          lastCoord = streamlines[loopB].yCoords[streamlines[loopB].totalPoints-1];
          break;
        case kdirZ: 
          lastCoord = streamlines[loopB].zCoords[streamlines[loopB].totalPoints-1];
          break;
      }      
      if ((fabs(lastCoord-maxCoord)<1.0e-2*refLength)&&(streamlines[loopB].totalPoints>=loopA)){
        // Update Counter
        count++;
        // Get Local Coords
        switch(dir){
          case kdirX: 
            xDiff = streamlines[loopB].yCoords[loopA]-streamlines[loopB].yCoords[0];
            yDiff = streamlines[loopB].zCoords[loopA]-streamlines[loopB].zCoords[0];
            break;
          case kdirY: 
            xDiff = streamlines[loopB].zCoords[loopA]-streamlines[loopB].zCoords[0];
            yDiff = streamlines[loopB].xCoords[loopA]-streamlines[loopB].xCoords[0];            
            break;
          case kdirZ: 
            xDiff = streamlines[loopB].xCoords[loopA]-streamlines[loopB].xCoords[0];
            yDiff = streamlines[loopB].yCoords[loopA]-streamlines[loopB].yCoords[0];            
            break;
        }              
        // Get Distance In plane
        dist = sqrt(xDiff*xDiff+yDiff*yDiff);
        // Add Contribution
        currentAv = currentAv + dist;
        currentSqrAv = currentSqrAv + dist*dist;
      }
    }
    // Normalize
    // Av
    if (count>0) currentAv = (currentAv/count);
    else currentAv = 0.0;
    // Sqr Av
    if (count>0) currentSqrAv = (currentSqrAv/count);
    else currentSqrAv = 0.0;
    value = currentSqrAv - currentAv*currentAv;
    if (value>kMathZero) crossDeviations[loopA] = sqrt(value);
    else crossDeviations[loopA] = 0.0;
  }
  // Eval The Number of Finishing Streamlines
  count = 0;
  for(int loopB=0;loopB<totalSL;loopB++){
    switch(dir){
      case kdirX: 
        lastCoord = streamlines[loopB].xCoords[streamlines[loopB].totalPoints];
        break;
      case kdirY: 
        lastCoord = streamlines[loopB].yCoords[streamlines[loopB].totalPoints];
        break;
      case kdirZ: 
        lastCoord = streamlines[loopB].zCoords[streamlines[loopB].totalPoints];
        break;
    }              
    if (fabs(lastCoord - maxCoord)<1.0e-2*refLength) {
      count++;
    }
  }
  // Message
  printf("Cross Diffusion Statistics Based On %d samples.\n",count);
}

// Eval Transverse Diffusion Statistics
void MRIScan::EvalSLTransverseDiffusionWithSpace(int totalSL, std::vector<MRIStreamline*> &streamlines,
                                                 MRIDirection dir, double minCoord, double maxCoord,
                                                 int totalSteps, std::vector<double> &space, std::vector<double> &crossDeviations){
  double currentAv = 0.0;
  double currentSqrAv = 0.0;
  double currentCoord = 0.0;
  double lastCoord = 0.0;
  double xDiff = 0.0;
  double yDiff = 0.0;
  double value = 0.0;
  int count = 0;
  bool success = false;
  double dist = 0.0;
  double* sLInt = new double[2];
  // Initialize
  double refLength = maxCoord - minCoord;
  double deltaZ = (double)refLength/(double)totalSteps;  
  // Allocate and Initialize
  space.reserve(totalSteps);
  crossDeviations.reserve(totalSteps);
  space[0] = minCoord;
  for(int loopA=1;loopA<totalSteps;loopA++){
    space[loopA] = space[loopA-1] + deltaZ;
  }
  // Loop On Time Steps
  for(int loopA=0;loopA<totalSteps;loopA++){
    currentAv = 0.0;
    currentSqrAv = 0.0;
    count = 0;
    currentCoord = space[loopA];
    for(int loopB=0;loopB<totalSL;loopB++){
      switch(dir){
        case kdirX: 
          lastCoord = streamlines[loopB]->xCoords[streamlines[loopB]->totalPoints-1];
          break;
        case kdirY: 
          lastCoord = streamlines[loopB]->yCoords[streamlines[loopB]->totalPoints-1];
          break;
        case kdirZ: 
          lastCoord = streamlines[loopB]->zCoords[streamlines[loopB]->totalPoints-1];
          break;
      }              
      if ((fabs(lastCoord-maxCoord)<1.0e-2*refLength)&&(streamlines[loopB]->totalPoints>=loopA)){
        success = streamlines[loopB]->EvalSLIntersection(dir,currentCoord,sLInt);
        if (success){
          // Update Counter
          count++;
          // Get Distance In plane          
          switch(dir){
            case kdirX: 
              xDiff = sLInt[0]-streamlines[loopB]->yCoords[0];
              yDiff = sLInt[1]-streamlines[loopB]->zCoords[0];
              break;
            case kdirY: 
              xDiff = sLInt[0]-streamlines[loopB]->zCoords[0];
              yDiff = sLInt[1]-streamlines[loopB]->xCoords[0];
              break;
            case kdirZ: 
              xDiff = sLInt[0]-streamlines[loopB]->xCoords[0];
              yDiff = sLInt[1]-streamlines[loopB]->yCoords[0];
              break;
          }              
          dist = sqrt(xDiff*xDiff+yDiff*yDiff);
          // Aggiungi
          currentAv = currentAv + dist;
          currentSqrAv = currentSqrAv + dist*dist;
        }
      }
    }
    // Normalize
    // Av
    if (count>0){
      currentAv = (currentAv/count);
    }else{
      currentAv = 0.0;
    }
    // Sqr Av
    if (count>0){
      currentSqrAv = (currentSqrAv/count);
    }else{
      currentSqrAv = 0.0;
    }
    value = currentSqrAv - currentAv*currentAv;
    if (value>kMathZero){
      crossDeviations[loopA] = sqrt(value);
    }else{
      crossDeviations[loopA] = 0.0;
    }
  }
  // Eval The Number of Finishing Streamlines
  count = 0;
  for(int loopB=0;loopB<totalSL;loopB++){
    lastCoord = streamlines[loopB]->zCoords[streamlines[loopB]->totalPoints-1];
    if (fabs(lastCoord-maxCoord)<1.0e-2*refLength){
      count++;
    }
  }
  // Message
  printf("Cross Diffusion Statistics Based On %d samples.\n",count);
}

// EVAL STREAMLINES STATISTICS
void MRIScan::EvalStreamLineStatistics(MRIDirection dir, MRIStreamlineOptions &options, std::vector<MRIStreamline*> &streamlines){
  double maxCoord = 0.0;
  double minCoord = 0.0;
  std::vector<double> sliceCenter;
  std::vector<double> time;
  std::vector<double> sliceNormArrivals;
  std::vector<double> crossDeviations;
  // Set Parameters
  int totalSL = options.gridTotals[0]*options.gridTotals[1];
  // Check The Direction
  switch(dir){
    case kdirX:
      maxCoord = domainSizeMax[0];
      minCoord = domainSizeMin[0];
      break;
    case kdirY:
      maxCoord = domainSizeMax[1];
      minCoord = domainSizeMin[1];
      break;    
    case kdirZ:
      maxCoord = domainSizeMax[2];
      minCoord = domainSizeMin[2];
      break;
  }
  // SET PARAMETERS
  int totalTimeSteps = (MRIUtils::round(options.totalT/options.deltaT)+1);
  int totalSlices = 50;
  // Arrival Point Distribution
  EvalSLArrivalPointDistribution(totalSL,streamlines,dir,minCoord,maxCoord,totalSlices,sliceCenter,sliceNormArrivals);
  // Write To File
  MRIUtils::WriteGraphToFile("SLArrivals.dat",totalSlices,sliceCenter,sliceNormArrivals);
  // Transverse Diffusion
  EvalSLTransverseDiffusionWithSpace(totalSL,streamlines,dir,minCoord,maxCoord,totalTimeSteps,time,crossDeviations);
  // Write To File
   MRIUtils::WriteGraphToFile("CrossDeviation.dat",totalTimeSteps,time,crossDeviations);
}

// Print StreamLines To File
void MRIStreamline::AppendToFile(int index, std::string fileName){
    // Open Output File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"a");
  // Write vertices
  for(int loopA=0;loopA<totalPoints;loopA++){
    fprintf(outFile,"%e,%e,%e,%d\n",xCoords[loopA],yCoords[loopA],zCoords[loopA],index);
  }
  // Write line definition
  //fprintf(outFile,"l ");
  //for(int loopA=0;loopA<totalPoints;loopA++){
  //  fprintf(outFile,"%d ",-totalPoints+loopA);
  //}
  //fprintf(outFile,"\n");

  // Close Output file
  fclose(outFile);
}
