# include "mriIO.h"

using namespace std;

// =======================
// READ SCAN FROM PLT FILE
// =======================
void readPLTData(string PltFileName, MRITopology* topo, MRIScan* scan){
  // Init Line Count
  int lineCount = 0;
  topology->totalCells = 0;

  // Initialize Plt Option Record
  PLTOptionRecord pltOptions;

  // Init Domain Limits
  topology->domainSizeMin[0] =  std::numeric_limits<double>::max();
  topology->domainSizeMin[1] =  std::numeric_limits<double>::max();
  topology->domainSizeMin[2] =  std::numeric_limits<double>::max();
  topology->domainSizeMax[0] = -std::numeric_limits<double>::max();
  topology->domainSizeMax[1] = -std::numeric_limits<double>::max();
  topology->domainSizeMax[2] = -std::numeric_limits<double>::max();

  // Assign File
  std::ifstream PltFile;
  writeSchMessage(std::string("Open File: ") + PltFileName + std::string("\n"));
  PltFile.open(PltFileName.c_str());

  // Init
  int TotalXCoords = 0;
  int TotalYCoords = 0;
  int TotalZCoords = 0;
  maxVelModule = 0.0;

  // Read The Number Of Lines
  std::vector<std::string> tokenizedString;
  bool foundheader = false;
  bool areAllFloats = false;
  int headerCount = 0;
  writeSchMessage(std::string("Computing input file size...\n"));
  int totalLinesInFile = 0;
  std::string Buffer;
  while (std::getline(PltFile,Buffer)){
    if(!foundheader){
      boost::trim(Buffer);
      boost::split(tokenizedString, Buffer, boost::is_any_of("= ,"), boost::token_compress_on);
      areAllFloats = true;
      assignPLTOptions(tokenizedString, pltOptions);
      for(size_t loopA=0;loopA<tokenizedString.size();loopA++){
        areAllFloats = (areAllFloats && (MRIUtils::isFloat(tokenizedString[loopA])));
      }
      foundheader = areAllFloats;
      headerCount++;
    }
    // Increase cell and line number
    totalLinesInFile++;
  }
  
  // Done: Computing Input File Size
  writeSchMessage(string("Header Size: " + MRIUtils::intToStr(headerCount) + "\n"));
  writeSchMessage(string("Total Lines: " + MRIUtils::intToStr(totalLinesInFile) + "\n"));
  writeSchMessage(std::string("Done.\n"));

  // Reset File
  PltFile.clear();
  PltFile.seekg(0, std::ios::beg);
  MRICell myCellPoint;

  // Skip Comments
  std::string* PltFileHeader = new std::string[headerCount];
  for(int loopA=0;loopA<headerCount-1;loopA++){
    std::getline(PltFile,Buffer);
    PltFileHeader[loopA] = Buffer;
    lineCount++;
  }

  // Initialize Local Variables   
  double LocalXCoord = 0.0;
  double LocalYCoord = 0.0;
  double LocalZCoord = 0.0;
  double LocalConc = 0.0;
  double LocalXVel = 0.0;
  double LocalYVel = 0.0;
  double LocalZVel = 0.0;
  double CurrentModule = 0.0;
  bool Continue = false;
  std::string outString = "";
  // Vector with X,Y,Z Coords
  MRIDoubleVec XCoords;
  MRIDoubleVec YCoords;
  MRIDoubleVec ZCoords;
  
  // Read All Lines
  int LocalCount = 0;
  int precentProgress = 0;
  int percentCounted = 0;
  int valueCounter = 0;
  int neededValues = 7;
  double* LocalVal = new double[neededValues];
  double* TempVal = new double[neededValues];
  vector<string> ResultArray;
  
  // Reading Input File Message
  writeSchMessage(std::string("Reading input file...\n"));
  
  while (std::getline(PltFile,Buffer)){
    // Read Line
    lineCount++;
    precentProgress = (int)(((double)lineCount/(double)totalLinesInFile)*100);
    if (((precentProgress % 10) == 0)&&((precentProgress / 10) != percentCounted)){
      percentCounted = (precentProgress / 10);
      writeSchMessage(std::string("Reading..."+MRIUtils::intToStr(precentProgress)+"\n"));
    }

    // Tokenize Line
    ResultArray = MRIUtils::extractSubStringFromBufferMS(Buffer);
    
    // Store Local Structure
    try{
      // Set Continue
      Continue = true;
      // Check Ratio between ResultArray.size, valueCounter, neededValues
      if((int)ResultArray.size()+valueCounter < neededValues){
        // Read the whole Result Array
        for(size_t loopA=0;loopA<ResultArray.size();loopA++){
          LocalVal[loopA] = atof(ResultArray[loopA].c_str());
        }
        Continue = false;
      }else{
        // Read part of the result array
        for(int loopA=0;loopA<neededValues-valueCounter;loopA++){
          LocalVal[loopA] = atof(ResultArray[loopA].c_str());
        }
        // Put the rest in temporary array
        for(size_t loopA=0;loopA<ResultArray.size()-(neededValues-valueCounter);loopA++){
          TempVal[loopA] = atof(ResultArray[loopA].c_str());
        }
        Continue = true;
        // Coords
        LocalXCoord = LocalVal[0];
        LocalYCoord = LocalVal[1];
        LocalZCoord = LocalVal[2];
        // Concentration
        LocalConc = LocalVal[3];
        // Velocity
        LocalXVel = LocalVal[4];
        LocalYVel = LocalVal[5];
        LocalZVel = LocalVal[6];
      }
      Continue = true;
      // Coords
      LocalXCoord = LocalVal[0];
      LocalYCoord = LocalVal[1];
      LocalZCoord = LocalVal[2];
      // Concentration
      LocalConc = LocalVal[3];
      // Velocity
      LocalXVel = LocalVal[4];
      LocalYVel = LocalVal[5];
      LocalZVel = LocalVal[6];
      // Update valueCounter
      valueCounter = ((ResultArray.size() + valueCounter) % neededValues);
      // Check Module
      CurrentModule = sqrt((LocalXVel*LocalXVel)+(LocalYVel*LocalYVel)+(LocalZVel*LocalZVel));
      if (CurrentModule>1000.0){
        throw 20;
      }
    }catch (...){
      //Set Continue
      Continue = false;
      std::string outString = "WARNING[*] Error Reading Line: "+MRIUtils::intToStr(lineCount)+"; Line Skipped.\n";
      printf("%s",outString.c_str());
    }
    if (Continue){
      // Update Limits
      // Min
      if (LocalXCoord<topology->domainSizeMin[0]) topology->domainSizeMin[0] = LocalXCoord;
      if (LocalYCoord<topology->domainSizeMin[1]) topology->domainSizeMin[1] = LocalYCoord;
      if (LocalZCoord<topology->domainSizeMin[2]) topology->domainSizeMin[2] = LocalZCoord;
      // Max
      if (LocalXCoord>topology->domainSizeMax[0]) topology->domainSizeMax[0] = LocalXCoord;
      if (LocalYCoord>topology->domainSizeMax[1]) topology->domainSizeMax[1] = LocalYCoord;
      if (LocalZCoord>topology->domainSizeMax[2]) topology->domainSizeMax[2] = LocalZCoord;

      // Update Max Speeds
      if (CurrentModule>maxVelModule) {
        maxVelModule = CurrentModule;
      }

      // Store Node Coords To Find Grid Size
      MRIUtils::insertInList(LocalXCoord,XCoords);
      MRIUtils::insertInList(LocalYCoord,YCoords);
      MRIUtils::insertInList(LocalZCoord,ZCoords);    

      // Store Velocity/Concentrations
      LocalCount++;
    
      // Conc
      myCellPoint.concentration = LocalConc;
      // Velocity
      myCellPoint.velocity[0] = LocalXVel;
      myCellPoint.velocity[1] = LocalYVel;
      myCellPoint.velocity[2] = LocalZVel;
    
      // Add to Vector
      cells.push_back(myCellPoint);

      // Set Continue
      Continue = true;
    }
  }
  delete [] LocalVal;
  delete [] TempVal;

  // Set The Effective Number Of Data Read
  topology->totalCells = LocalCount;

  // Store Total Cells
  topology->cellTotals[0] = XCoords.size();
  topology->cellTotals[1] = YCoords.size();
  topology->cellTotals[2] = ZCoords.size();

  if(topology->totalCells != TotalXCoords * TotalYCoords * TotalZCoords){
    writeSchMessage(std::string("Total number of cells in X: " + MRIUtils::intToStr(TotalXCoords) + "\n"));
    writeSchMessage(std::string("Total number of cells in Y: " + MRIUtils::intToStr(TotalYCoords) + "\n"));
    writeSchMessage(std::string("Total number of cells in Z: " + MRIUtils::intToStr(TotalZCoords) + "\n"));
    writeSchMessage(std::string("Total number of cells: " + MRIUtils::intToStr(topology->totalCells) + "\n"));
    throw MRIException("ERROR: Total Number of Cells does not match!\n");
  }

  // Complete To Full Grid: Set To Zero
  topology->totalCells = TotalXCoords * TotalYCoords * TotalZCoords;

  // Set a Zero mtCellPoint
  myCellPoint.concentration = 0.0;
  myCellPoint.velocity[0] = 0.0;
  myCellPoint.velocity[1] = 0.0;
  myCellPoint.velocity[2] = 0.0;

  // Resize CellPoints
  cells.resize(topology->totalCells,myCellPoint);

  // Resize CellLenghts: UNIFORM CASE
  topology->cellLengths.resize(3);
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    topology->cellLengths[loopA].resize(topology->cellTotals[loopA]);
    if(topology->cellTotals[loopA] == 1){
      topology->cellLengths[loopA][0] = 1.0;
    }else{
      for(int loopB=0;loopB<topology->cellTotals[loopA];loopB++){
        topology->cellLengths[loopA][loopB] = fabs(topology->domainSizeMax[loopA]-topology->domainSizeMin[loopA])/(topology->cellTotals[loopA]-1);
      }
    }
  }

  // Finished Reading File
  writeSchMessage(std::string("File reading completed.\n"));

  // Close File
  PltFile.close();

  // REORDER CELLS
  if (DoReorderCells){
    reorderScan();
  }

  // WRITE STATISTICS
  std::string CurrentStats = WriteStatistics();
  writeSchMessage(CurrentStats);

}
