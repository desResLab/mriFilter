# include "mriScan.h"

using namespace std;

MRIScan::MRIScan(double currentTime){
  scanTime = currentTime;
}

// ==============================================
// READS TOKENIZED STRING AND ASSIGNS PLT OPTIONS
// ==============================================
void assignVTKOptions(int lineNum, std::vector<std::string> tokens, vtkStructuredPointsOptionRecord &vtkOptions){
  for(size_t loopA=0;loopA<tokens.size();loopA++){
    if(boost::to_upper_copy(tokens[loopA]) == "ASCII"){
      vtkOptions.isASCII = true;
      vtkOptions.isDefined[0] = true;
    }else if(boost::to_upper_copy(tokens[loopA]).find("STRUCTURED") != string::npos){
      vtkOptions.isValidDataset = true;
      vtkOptions.isDefined[1] = true;
    }else if(boost::to_upper_copy(tokens[loopA]) == "DIMENSIONS"){
      vtkOptions.dimensions[0] = atoi(tokens[loopA+1].c_str());
      vtkOptions.dimensions[1] = atoi(tokens[loopA+2].c_str());
      vtkOptions.dimensions[2] = atoi(tokens[loopA+3].c_str());
      vtkOptions.isDefined[2] = true;
    }else if(boost::to_upper_copy(tokens[loopA]) == "ORIGIN"){
      vtkOptions.origin[0] = atof(tokens[loopA+1].c_str());
      vtkOptions.origin[1] = atof(tokens[loopA+2].c_str());
      vtkOptions.origin[2] = atof(tokens[loopA+3].c_str());
      vtkOptions.isDefined[3] = true;
    }else if(boost::to_upper_copy(tokens[loopA]) == "SPACING"){
      vtkOptions.spacing[0] = atof(tokens[loopA+1].c_str());
      vtkOptions.spacing[1] = atof(tokens[loopA+2].c_str());
      vtkOptions.spacing[2] = atof(tokens[loopA+3].c_str());
      vtkOptions.isDefined[4] = true;
    }else if(boost::to_upper_copy(tokens[loopA]) == "SCALARS"){
      // Add Starting Point
      vtkOptions.dataBlockStart.push_back(lineNum);
      // Add type
      vtkOptions.dataBlockType.push_back(0);
      if((boost::to_upper_copy(tokens[loopA + 1]) == "CONCENTRATION")){
        vtkOptions.dataBlockRead.push_back(true);
      }else{
        vtkOptions.dataBlockRead.push_back(false);
      }
    }else if(boost::to_upper_copy(tokens[loopA]) == "VECTORS"){
      // Add Starting Point
      vtkOptions.dataBlockStart.push_back(lineNum);
      // Add type
      vtkOptions.dataBlockType.push_back(1);
      if((boost::to_upper_copy(tokens[loopA + 1]) == "VELOCITY")){
        vtkOptions.dataBlockRead.push_back(true);
      }else{
        vtkOptions.dataBlockRead.push_back(false);
      }
    }
  }
}

// ================
// COPY CONSTRUCTOR
// ================
MRIScan::MRIScan(const MRIScan& copyScan){
  // Assign Scan Time
  scanTime = copyScan.scanTime;
  // Resent Pressure Gradient and Relative Pressure
  hasPressureGradient = copyScan.hasPressureGradient;
  hasRelativePressure = copyScan.hasRelativePressure;
  hasReynoldsStress = copyScan.hasReynoldsStress;
}

// Print the File List Log
void PrintFileListLog(int totalFiles,std::string* fileNames){
  // Open Output File
  FILE* outFile;
  outFile = fopen("FileList.log","w");
  // Write Header
  for(int loopA=0;loopA<totalFiles;loopA++){
    fprintf(outFile,"%d %s\n",loopA+1,fileNames[loopA].c_str());
  }
  // Close Output file
  fclose(outFile);
}

// ========================
// FILL TECPLOT FILE HEADER
// ========================
void MRIScan::fillPLTHeader(std::vector<std::string> &pltHeader, bool isFirstFile){
  // Clear Vector
  pltHeader.clear();
  if (isFirstFile){
    pltHeader.push_back("TITLE = \"smpFilterOutput\"");
    pltHeader.push_back("VARIABLES = \"X/D\"");
    pltHeader.push_back("\"Y/D\"");
    pltHeader.push_back("\"Z/D\"");
    pltHeader.push_back("\"Conc%\"");
    pltHeader.push_back("\"Vx/Ubulk\"");
    pltHeader.push_back("\"Vy/Ubulk\"");
    pltHeader.push_back("\"Vz/Ubulk\"");
    // Print all outputs
    for(size_t loopA=0;loopA<outputs.size();loopA++){
      // Add output Name
      if(outputs[loopA].totComponents == 1){
        pltHeader.push_back("\"" + outputs[loopA].name + "\"");
      }else if(outputs[loopA].totComponents == 2){
        pltHeader.push_back("\"" + outputs[loopA].name + "x\"");
        pltHeader.push_back("\"" + outputs[loopA].name + "y\"");
      }else if(outputs[loopA].totComponents == 3){
        pltHeader.push_back("\"" + outputs[loopA].name + "x\"");
        pltHeader.push_back("\"" + outputs[loopA].name + "y\"");
        pltHeader.push_back("\"" + outputs[loopA].name + "z\"");
      }else{
        throw MRIException("Too many components for TECPLOT result.\n");
      }
    }
  }
  pltHeader.push_back("ZONE T=\"SubZone\"");
  pltHeader.push_back(" STRANDID=0, SOLUTIONTIME="+MRIUtils::floatToStr(scanTime));
  pltHeader.push_back(" I=35, J=113, K=155, ZONETYPE=Ordered");
  pltHeader.push_back(" DATAPACKING=POINT");
  string singleString = string(" DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE ");
  for(size_t loopA=0;loopA<outputs.size();loopA++){
    for(int loopB=0;loopB<outputs[loopA].totComponents;loopB++){
      singleString += "SINGLE ";
    }
  }
  singleString += ")";
  pltHeader.push_back(singleString);
}

// returns file size in bytes or -1 if not found.
std::ifstream::pos_type GetFileSize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::in | std::ifstream::binary);
    in.seekg(0, std::ifstream::end);
    return in.tellg(); 
}

// Write IO Log
void WriteIOLog(std::string LogFileName, std::string MsgsString)
{
  // Open Output File
	FILE* outFile;
	outFile = fopen(LogFileName.c_str(),"a");
	// Write Header
  fprintf(outFile,"%s\n",MsgsString.c_str());
	// Close Output file
	fclose(outFile);	
}

// =============
// REORDER CELLS
// =============
void MRIScan::reorderScan(){
  
  // Determine The Direct and Inverse Permutations
  WriteSchMessage(std::string("Computing Permutation..."));
  std::vector<int> DirectPerm;
  getGlobalPermutation(DirectPerm);
  WriteSchMessage(std::string("Done.\n"));
  
  // Reorder Cells
  WriteSchMessage(std::string("Reordering Cells..."));
  reorderCells(DirectPerm);
  WriteSchMessage(std::string("Done.\n"));
}

// ====================================
// READS HEADER AND ASSIGNS PLT OPTIONS
// ====================================
void assignPLTOptions(std::vector<std::string> tokens, PLTOptionRecord &pltOptions){
  for(size_t loopA=0;loopA<tokens.size();loopA++){
    if(boost::to_upper_copy(tokens[loopA]) == "I"){
      pltOptions.type = pltUNIFORM;
      pltOptions.i = atoi(tokens[loopA+1].c_str());
    }else if(boost::to_upper_copy(tokens[loopA]) == "J"){
      pltOptions.type = pltUNIFORM;
      pltOptions.j = atoi(tokens[loopA+1].c_str());
    }else if(boost::to_upper_copy(tokens[loopA]) == "K"){
      pltOptions.type = pltUNIFORM;
      pltOptions.k = atoi(tokens[loopA+1].c_str());
    }else if(boost::to_upper_copy(tokens[loopA]) == "N"){
      pltOptions.type = pltSTRUCTURED;
      pltOptions.N = atoi(tokens[loopA+1].c_str());
    }else if(boost::to_upper_copy(tokens[loopA]) == "E"){
      pltOptions.type = pltSTRUCTURED;
      pltOptions.E = atoi(tokens[loopA+1].c_str());
    }else if(boost::to_upper_copy(tokens[loopA]) == "FEBLOCK"){
      pltOptions.type = pltSTRUCTURED;
    }else if(boost::to_upper_copy(tokens[loopA]) == "DATAPACKING"){
      if(boost::to_upper_copy(tokens[loopA+1]) != "POINT"){
        throw MRIException("Error: invalid DATAPACKING format");
      }
    }
  }
}

// =======================
// READ SCAN FROM PLT FILE
// =======================
void MRIScan::readPLTFile(std::string PltFileName, bool DoReorderCells){
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
  WriteSchMessage(std::string("Open File: ") + PltFileName + std::string("\n"));
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
  WriteSchMessage(std::string("Computing input file size...\n"));
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
  WriteSchMessage(string("Header Size: " + MRIUtils::intToStr(headerCount) + "\n"));
  WriteSchMessage(string("Total Lines: " + MRIUtils::intToStr(totalLinesInFile) + "\n"));
  WriteSchMessage(std::string("Done.\n"));

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
  WriteSchMessage(std::string("Reading input file...\n"));
  
  while (std::getline(PltFile,Buffer)){
    // Read Line
    lineCount++;
    precentProgress = (int)(((double)lineCount/(double)totalLinesInFile)*100);
    if (((precentProgress % 10) == 0)&&((precentProgress / 10) != percentCounted)){
      percentCounted = (precentProgress / 10);
      WriteSchMessage(std::string("Reading..."+MRIUtils::intToStr(precentProgress)+"\n"));
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
    WriteSchMessage(std::string("Total number of cells in X: " + MRIUtils::intToStr(TotalXCoords) + "\n"));
    WriteSchMessage(std::string("Total number of cells in Y: " + MRIUtils::intToStr(TotalYCoords) + "\n"));
    WriteSchMessage(std::string("Total number of cells in Z: " + MRIUtils::intToStr(TotalZCoords) + "\n"));
    WriteSchMessage(std::string("Total number of cells: " + MRIUtils::intToStr(topology->totalCells) + "\n"));
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
  WriteSchMessage(std::string("File reading completed.\n"));

  // Close File
  PltFile.close();

  // REORDER CELLS
  if (DoReorderCells){
    reorderScan();
  }

  // WRITE STATISTICS
  std::string CurrentStats = WriteStatistics();
  WriteSchMessage(CurrentStats);

}

// ================
// EXPORT TO LSDYNA
// ================
void MRIScan::exportToLSDYNA(std::string LSFileName, double scale){
  // Export to LS-DYNA
  printf("Exporting File to LSDyna...");
  FILE* LSFile;
  LSFile = fopen(LSFileName.c_str(),"w");
  // Write Header
  fprintf(LSFile,"*KEYWORD\n");
  fprintf(LSFile,"*NODE\n");

  // Create All Nodes
  double CurrentXCoord,CurrentYCoord,CurrentZCoord;
  int TotalNodes = 0;
	double CurrentModule;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    CurrentModule = sqrt((cells[loopA].velocity[0]*cells[loopA].velocity[0])+
                         (cells[loopA].velocity[1]*cells[loopA].velocity[1])+
                         (cells[loopA].velocity[2]*cells[loopA].velocity[2]));
    // Export Nodes Based On Original Velocity
    if (CurrentModule>kMathZero){
      // Node 1
      TotalNodes++;
      CurrentXCoord = topology->cellLocations[loopA][0] - cells[loopA].velocity[0]*scale;
      CurrentYCoord = topology->cellLocations[loopA][1] - cells[loopA].velocity[1]*scale;
      CurrentZCoord = topology->cellLocations[loopA][2] - cells[loopA].velocity[2]*scale;
      fprintf(LSFile,"%d,%e,%e,%e,0,0\n",TotalNodes,CurrentXCoord,CurrentYCoord,CurrentZCoord);
      // Node 2
      TotalNodes++;
      CurrentXCoord = topology->cellLocations[loopA][0] + cells[loopA].velocity[0]*scale;
      CurrentYCoord = topology->cellLocations[loopA][1] + cells[loopA].velocity[1]*scale;
      CurrentZCoord = topology->cellLocations[loopA][2] + cells[loopA].velocity[2]*scale;
      fprintf(LSFile,"%d,%e,%e,%e,0,0\n",TotalNodes,CurrentXCoord,CurrentYCoord,CurrentZCoord);
    }
  }
  // WRITE ELEMENTS
  fprintf(LSFile,"*ELEMENT_BEAM\n");
  int Count = 0;
	int Node1,Node2;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    CurrentModule = sqrt((cells[loopA].velocity[0]*cells[loopA].velocity[0])+
                         (cells[loopA].velocity[1]*cells[loopA].velocity[1])+
                         (cells[loopA].velocity[2]*cells[loopA].velocity[2]));
    if (CurrentModule>kMathZero){
      Count++;
      Node1 = (Count-1)*2+1;
      Node2 = Count*2;
      fprintf(LSFile,"%d,1,%d,%d\n",Count,Node1,Node2);
    }
  }
  // Close File
  fclose(LSFile);
  printf("Done\n");	
};

// =============
// EXPORT TO CSV
// =============
void MRIScan::exportToCSV(std::string FileName){
  printf("Exporting to CSV...");
  std::ofstream OutFile;
  OutFile.open(FileName.c_str());
  // Loop On Cells
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    OutFile << topology->cellLocations[loopA][0] << 
               topology->cellLocations[loopA][1] << 
               topology->cellLocations[loopA][2] <<
               cells[loopA].concentration <<
			         cells[loopA].velocity[0] << 
               cells[loopA].velocity[1] << 
               cells[loopA].velocity[2];
  }
  OutFile.close();
  printf("Done\n");
}

// ============================
// EXPORT TO TECPLOT ASCII FILE
// ============================
void MRIScan::exportToTECPLOT(std::string FileName, bool isFirstFile){
  // Write Progress Message 
	WriteSchMessage(std::string("Exporting to TECPLOT..."));
  
  // Open Output File
	FILE* outFile;
  if (isFirstFile){
	  outFile = fopen(FileName.c_str(),"w");
  }else{
    outFile = fopen(FileName.c_str(),"a");
  }

  // Fill Header
  std::vector<std::string> PltFileHeader;
  fillPLTHeader(PltFileHeader,isFirstFile);

  // Write Header
  std::string LineString;
  std::string compString = "I";
  for(unsigned int loopA=0;loopA<PltFileHeader.size();loopA++){
    boost::trim(PltFileHeader[loopA]);
    LineString = PltFileHeader[loopA];
    if (LineString.substr(0,1) != compString){
      fprintf(outFile,"%s\n",PltFileHeader[loopA].c_str());
    }else{
      fprintf(outFile," I=%d, J=%d, K=%d, ZONETYPE=Ordered\n",topology->cellTotals[0],topology->cellTotals[1],topology->cellTotals[2]);
    }
  }
  // Loop On Cells
  for(int loopA=0;loopA<topology->totalCells;loopA++){

    // Write position, concentration and velocity
    fprintf(outFile,"%-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e ",
                      topology->cellLocations[loopA][0],
                      topology->cellLocations[loopA][1],
                      topology->cellLocations[loopA][2],
                      cells[loopA].concentration,
                      cells[loopA].velocity[0],
                      cells[loopA].velocity[1],
                      cells[loopA].velocity[2]);

    // Add result quantities
    for(size_t loopB=0;loopB<outputs.size();loopB++){
      for(int loopC=0;loopC<outputs[loopB].totComponents;loopC++){
        fprintf(outFile,"%-15.6e ",outputs[loopB].values[loopA*outputs[loopB].totComponents + loopC]);
      }
    }

    // New Line
    fprintf(outFile,"\n");
  }

  // Close Output file
  fclose(outFile);
  
  // Write Done Message
  WriteSchMessage(std::string("Done\n"));
};

// ========================
// Get Local Adjacent Plane
// ========================
void MRIScan::getLocalStarFaces(int StarNum, int CellsX, int CellsY, int &BottomFace, int &TopFace, int &LeftFace, int &RightFace)
{
  // Find Local Face Number
  BottomFace = (((int)(StarNum) / (int)(CellsX+1))-1)*(2*CellsX+1) + ((int)(StarNum) % (int)(CellsX+1)) + CellsX;
  TopFace =    (((int)(StarNum) / (int)(CellsX+1)))  *(2*CellsX+1) + ((int)(StarNum) % (int)(CellsX+1)) + CellsX;
  LeftFace =   (((int)(StarNum) / (int)(CellsX+1)))  *(2*CellsX+1) + ((int)(StarNum) % (int)(CellsX+1)) - 1;
  RightFace =  LeftFace + 1;

  // Set To Zero the Null Faces
  int iCoord = ((int)(StarNum) / (int)(CellsX+1));
  int jCoord = ((int)(StarNum) % (int)(CellsX+1));
  if (iCoord == 0) BottomFace = -1;
  if (iCoord == CellsY) TopFace = -1;
  if (jCoord == 0) LeftFace = -1;
  if (jCoord == CellsX) RightFace = -1;
};

// Flush Model To File
void MRIScan::flushToFile(std::string FileName){
  // Open Output File
	FILE* outFile;
	outFile = fopen(FileName.c_str(),"w");
	// Write Header
  fprintf(outFile,"%-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s\n","Cell Number","PosX","PosY","PosZ","Conc","VelX","VelY","VelZ");
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    fprintf(outFile,"%-15d %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n",loopA,
            topology->cellLocations[loopA][0],
            topology->cellLocations[loopA][1],
            topology->cellLocations[loopA][2],
            cells[loopA].concentration,
            cells[loopA].velocity[0],cells[loopA].velocity[1],cells[loopA].velocity[2]);
  }
	// Close Output file
	fclose(outFile);
}

// Eval The Central Cell for the Domain
int MRIScan::evalCentralCell(){
  return mapCoordsToIndex(topology->cellTotals[0]/2.0,topology->cellTotals[1]/2.0,topology->cellTotals[2]/2.0);
}

// ASSEMBLE ENCODING MATRIX
void MRIScan::assembleEncodingMatrix(int &totalRows, int &totalColumns, double** &Mat){
  // VAR
  int faceXPlus =  0;
  int faceXMinus = 0;
  int faceYPlus = 0;
  int faceYMinus = 0;
  int faceZPlus = 0;
  int faceZMinus = 0;
  int faceXColumn = 0;
  int faceYColumn = 0;
  int faceZColumn = 0;
  double faceXArea = 0.0;
  double faceYArea = 0.0;
  double faceZArea = 0.0;
  double Areas[3] = {0.0};

  
  // Values For Internal and Edges
  double edgeFactor = 0.5;
  double intFactor = 0.5;
  
  // FIND THE TOTAL NUMBER OF FACES
  totalRows = GetTotalFaces();
  totalColumns = 3*topology->totalCells;

  // ALLOCATE MATRIX
  int faceConn[totalRows];
  Mat = new double*[totalRows];
  for(int loopA=0;loopA<totalRows;loopA++){
    faceConn[loopA] = 0;
    Mat[loopA] = new double[totalColumns];
  }
  // INIT MATRIX
  for(int loopA=0;loopA<totalRows;loopA++){
    for(int loopB=0;loopB<totalColumns;loopB++){
      Mat[loopA][loopB] = 0.0;
    }
  }
  // FORM CONNECTIVITY MATRIX
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Eval Neighbours
    faceXPlus =  getAdjacentFace(loopA,kfacePlusX);
    faceXMinus = getAdjacentFace(loopA,kfaceMinusX);
    faceYPlus =  getAdjacentFace(loopA,kfacePlusY);
    faceYMinus = getAdjacentFace(loopA,kfaceMinusY);
    faceZPlus =  getAdjacentFace(loopA,kfacePlusZ);
    faceZMinus = getAdjacentFace(loopA,kfaceMinusZ);
    // Increment Counters
    faceConn[faceXPlus]++;
    faceConn[faceXMinus]++;
    faceConn[faceYPlus]++;
    faceConn[faceYMinus]++;
    faceConn[faceZPlus]++;
    faceConn[faceZMinus]++;
  }
  // LOOP THROUGH THE FACES
  for(int loopA=0;loopA<topology->totalCells;loopA++){

    // EVAL NEIGHBOURS
    faceXPlus =  getAdjacentFace(loopA,kfacePlusX);
    faceXMinus = getAdjacentFace(loopA,kfaceMinusX);
    faceYPlus =  getAdjacentFace(loopA,kfacePlusY);
    faceYMinus = getAdjacentFace(loopA,kfaceMinusY);
    faceZPlus =  getAdjacentFace(loopA,kfacePlusZ);
    faceZMinus = getAdjacentFace(loopA,kfaceMinusZ);

    // EVAL AREAS
    evalCellAreas(loopA,Areas);
    faceXArea = Areas[0];
    faceYArea = Areas[1];
    faceZArea = Areas[2];

    // Eval Column Number
    faceXColumn = loopA;
    faceYColumn = topology->totalCells + loopA;
    faceZColumn = 2*topology->totalCells + loopA;
    // Assembling Terms
    // X PLUS
    if(faceConn[faceXPlus] == 1){
      Mat[faceXPlus] [faceXColumn] += edgeFactor * faceXArea;
    }else if(faceConn[faceXPlus] == 2){
      Mat[faceXPlus] [faceXColumn] += intFactor * faceXArea;
    }else{
      throw MRIException("Invalid Face Connectivity");
    }
    // X MINUS
    if(faceConn[faceXMinus] == 1){
      Mat[faceXMinus][faceXColumn] += edgeFactor * faceXArea;
    }else if(faceConn[faceXMinus] == 2){
      Mat[faceXMinus][faceXColumn] += intFactor * faceXArea;
    }else{
      throw MRIException("Invalid Face Connectivity");
    }
    // Y PLUS
    if(faceConn[faceYPlus] == 1){
      Mat[faceYPlus] [faceYColumn] += edgeFactor * faceYArea;
    }else if(faceConn[faceYPlus] == 2){
      Mat[faceYPlus] [faceYColumn] += intFactor * faceYArea;
    }else{
      throw MRIException("Invalid Face Connectivity");
    }
    // Y MINUS
    if(faceConn[faceYMinus] == 1){
      Mat[faceYMinus][faceYColumn] += edgeFactor * faceYArea;
    }else if(faceConn[faceYMinus] == 2){
      Mat[faceYMinus][faceYColumn] += intFactor * faceYArea;
    }else{
      throw MRIException("Invalid Face Connectivity");
    }
    // Z PLUS
    if(faceConn[faceZPlus] == 1){
      Mat[faceZPlus] [faceZColumn] += edgeFactor * faceZArea;
    }else if(faceConn[faceZPlus] == 2){
      Mat[faceZPlus] [faceZColumn] += intFactor * faceZArea;
    }else{
      throw MRIException("Invalid Face Connectivity");
    }
    // Z MINUS
    if(faceConn[faceZMinus] == 1){
      Mat[faceZMinus][faceZColumn] += edgeFactor * faceZArea;
    }else if(faceConn[faceZMinus] == 2){
      Mat[faceZMinus][faceZColumn] += intFactor * faceZArea;
    }else{
      throw MRIException("Invalid Face Connectivity");
    }
  }
}

// Assemble Decoding Matrix
void MRIScan::assembleDecodingMatrix(int &totalRows, int &totalColumns, double** &Mat){
  // VAR
  int faceXPlus  = 0;
  int faceXMinus = 0;
  int faceYPlus  = 0;
  int faceYMinus = 0;
  int faceZPlus  = 0;
  int faceZMinus = 0;
  // Eval Column Number
  int faceXRow = 0;
  int faceYRow = 0;
  int faceZRow = 0;
  double Areas[3] = {0.0};
  double faceXArea = 0.0;
  double faceYArea = 0.0;
  double faceZArea = 0.0;
  
  // FIND THE TOTAL NUMBER OF FACES
  totalRows = 3*topology->totalCells;
  totalColumns = GetTotalFaces();
  
  // ALLOCATE MATRIX
  Mat = new double*[totalRows];
  for(int loopA=0;loopA<totalRows;loopA++){
    Mat[loopA] = new double[totalColumns];
  }
  // INIT MATRIX
  for(int loopA=0;loopA<totalRows;loopA++){
    for(int loopB=0;loopB<totalColumns;loopB++){
      Mat[loopA][loopB] = 0.0;
    }
  }
  // LOOP THROUGH THE FACES
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Get Areas
    evalCellAreas(loopA,Areas);
    faceXArea = Areas[0];
    faceYArea = Areas[1];
    faceZArea = Areas[2];

    // Eval Neighbours
    faceXPlus =  getAdjacentFace(loopA,kfacePlusX);
    faceXMinus = getAdjacentFace(loopA,kfaceMinusX);
    faceYPlus =  getAdjacentFace(loopA,kfacePlusY);
    faceYMinus = getAdjacentFace(loopA,kfaceMinusY);
    faceZPlus =  getAdjacentFace(loopA,kfacePlusZ);
    faceZMinus = getAdjacentFace(loopA,kfaceMinusZ);
    // Eval Column Number
    faceXRow = loopA;
    faceYRow = topology->totalCells + loopA;
    faceZRow = 2*topology->totalCells + loopA;
    // Assembling Terms
    Mat[faceXRow][faceXPlus]  += 0.5 * (1.0/faceXArea);
    Mat[faceXRow][faceXMinus] += 0.5 * (1.0/faceXArea);
    Mat[faceYRow][faceYPlus]  += 0.5 * (1.0/faceYArea);
    Mat[faceYRow][faceYMinus] += 0.5 * (1.0/faceYArea);
    Mat[faceZRow][faceZPlus]  += 0.5 * (1.0/faceZArea);
    Mat[faceZRow][faceZMinus] += 0.5 * (1.0/faceZArea);
  }
}

// Assign Random Component
void MRIScan::assignRandomComponent(const int direction,stdRndGenerator &generator){
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    switch (direction){
      case kdirX:
        cells[loopA].velocity[0] = generator();
        break;
      case kdirY:
        cells[loopA].velocity[1] = generator();
        break;
      case kdirZ:
        cells[loopA].velocity[2] = generator();
        break;
    }
  }
}

// Export Velocities To File As Row in the order X,Y,Z
void MRIScan::exportVelocitiesToFile(std::string fileName, bool append){
  
  // Open Output File
  FILE* outFile;
  
  // Check if Append
  if (append){
    outFile = fopen(fileName.c_str(),"a");
  }else{
    outFile = fopen(fileName.c_str(),"w");
  }
  // Write Header
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    for(int loopB=0;loopB<topology->totalCells;loopB++){
      fprintf(outFile,"%e ",cells[loopB].velocity[loopA]);  
    }
  }
  fprintf(outFile,"\n");  
  
  // Close Output file
  fclose(outFile);
}

// CROP SCAN
void MRIScan::crop(const MRIDoubleVec& limitBox){
  // Count The Number Of Cells Remaining
  int remainingCells = 0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    if (MRIUtils::isPointInsideBox(topology->cellLocations[loopA][0],
                                   topology->cellLocations[loopA][1],
                                   topology->cellLocations[loopA][2],
                                   limitBox)){
      remainingCells++;
    }
  }
  // Allocate New Cellpoints
  std::vector<MRICell> tempCellPoints;
  tempCellPoints.resize(remainingCells);
  maxVelModule = 0.0;
  int TotalXCoords = 0;
  int TotalYCoords = 0;
  int TotalZCoords = 0;
  std::vector<double> XCoords;
  std::vector<double> YCoords;
  std::vector<double> ZCoords;
  double currVelModule = 0.0;
  // Fill Temporary Cells
  int tempCount = 0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    if (MRIUtils::isPointInsideBox(topology->cellLocations[loopA][0],
                                   topology->cellLocations[loopA][1],
                                   topology->cellLocations[loopA][2],limitBox)){
      // Concentration
      tempCellPoints[tempCount].concentration = cells[loopA].concentration;
      // Velocity
      tempCellPoints[tempCount].velocity[0] = cells[loopA].velocity[0];
      tempCellPoints[tempCount].velocity[1] = cells[loopA].velocity[1];
      tempCellPoints[tempCount].velocity[2] = cells[loopA].velocity[2];
      // Get Velocity Module
      currVelModule = sqrt(tempCellPoints[tempCount].velocity[0]*tempCellPoints[tempCount].velocity[0]+
                           tempCellPoints[tempCount].velocity[1]*tempCellPoints[tempCount].velocity[1]+
                           tempCellPoints[tempCount].velocity[2]*tempCellPoints[tempCount].velocity[2]);
      // Store Maximum Velocity Module
      if (currVelModule>maxVelModule){
        maxVelModule = currVelModule;
      } 
      // Pressure Gradient
      tempCellPoints[tempCount].pressGrad[0] = cells[loopA].pressGrad[0];
      tempCellPoints[tempCount].pressGrad[1] = cells[loopA].pressGrad[1];
      tempCellPoints[tempCount].pressGrad[2] = cells[loopA].pressGrad[2];
      // Relative Pressures
      tempCellPoints[tempCount].relPressure = cells[loopA].relPressure;
      // Update Counter
      tempCount++;
    }
  }
  cells.resize(remainingCells);
  // Copy Back
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Concentration
    cells[loopA].concentration = tempCellPoints[loopA].concentration;
    // Velocity
    cells[loopA].velocity[0] = tempCellPoints[loopA].velocity[0];
    cells[loopA].velocity[1] = tempCellPoints[loopA].velocity[1];
    cells[loopA].velocity[2] = tempCellPoints[loopA].velocity[2];
    // Filtered Velocities
    cells[loopA].auxVector[0] = 0.0;
    cells[loopA].auxVector[1] = 0.0;
    cells[loopA].auxVector[2] = 0.0;
    // Pressure Gradients
    cells[loopA].pressGrad[0] = tempCellPoints[loopA].pressGrad[0];
    cells[loopA].pressGrad[1] = tempCellPoints[loopA].pressGrad[1];
    cells[loopA].pressGrad[2] = tempCellPoints[loopA].pressGrad[2];
    // Relative Pressure
    cells[loopA].relPressure = tempCellPoints[loopA].relPressure;
  }
}

// SCALE VELOCITIES
void MRIScan::scaleVelocities(double factor){
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    cells[loopA].velocity[0] *= factor;
    cells[loopA].velocity[1] *= factor;
    cells[loopA].velocity[2] *= factor;
  }
  maxVelModule *= factor;
}

// ===========================
// CHECK IF SPACING IS UNIFORM
// ===========================
bool MRIScan::hasUniformSpacing(){
  double lengthX = topology->cellLengths[0][0];
  double lengthY = topology->cellLengths[1][0];
  double lengthZ = topology->cellLengths[2][0];
  int count = 1;
  bool result = true;
  while(result && ((size_t)count<topology->cellLengths.size())){
    result = result && (fabs(topology->cellLengths[0][count] - lengthX) < kMathZero) &&
                       (fabs(topology->cellLengths[1][count] - lengthY) < kMathZero) &&
                       (fabs(topology->cellLengths[2][count] - lengthZ) < kMathZero);
    count++;
  }
  return result;
}

// =================
// Write to VTK File
// =================
void MRIScan::exportToVTK(std::string fileName, MRIThresholdCriteria* threshold){

  // Declare
  bool printAux = true;
  double currXCoord = 0.0;
  double currYCoord = 0.0;
  double currZCoord = 0.0;

  // Open Output File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  // Write Header
  fprintf(outFile,"# vtk DataFile Version 3.0\n");
  fprintf(outFile,"Grid Point Model\n");
  fprintf(outFile,"ASCII\n");

  if(hasUniformSpacing()){
      printf("Dataset type STRUCTURED_POINTS.\n");
      // Write Data Set
      fprintf(outFile,"DATASET STRUCTURED_POINTS\n");
      fprintf(outFile,"DIMENSIONS %d %d %d\n",topology->cellTotals[0],topology->cellTotals[1],topology->cellTotals[2]);
      fprintf(outFile,"SPACING %e %e %e\n",topology->cellLengths[0][0],topology->cellLengths[1][0],topology->cellLengths[2][0]);
      fprintf(outFile,"ORIGIN %e %e %e\n",topology->domainSizeMin[0],topology->domainSizeMin[1],topology->domainSizeMin[2]);
  }else{
    printf("Dataset type: RECTILINEAR_GRID.\n");
    // Write Data Set
    fprintf(outFile,"DATASET RECTILINEAR_GRID\n");
    fprintf(outFile,"DIMENSIONS %d %d %d\n",topology->cellTotals[0],topology->cellTotals[1],topology->cellTotals[2]);

    // Export X Coordinates
    fprintf(outFile,"X_COORDINATES %d double\n",(int)topology->cellLengths[0].size());
    currXCoord = topology->domainSizeMin[0];
    for(size_t loopA=1;loopA<topology->cellLengths[0].size();loopA++){
      fprintf(outFile,"%e ",currXCoord);
      currXCoord += 0.5*(topology->cellLengths[0][loopA-1] + topology->cellLengths[0][loopA]);
    }
    fprintf(outFile,"%e\n",currXCoord);

    // Export Y Coordinates
    fprintf(outFile,"Y_COORDINATES %d double\n",(int)topology->cellLengths[1].size());
    currYCoord = topology->domainSizeMin[1];
    for(size_t loopA=1;loopA<topology->cellLengths[1].size();loopA++){
      fprintf(outFile,"%e ",currYCoord);
      currYCoord += 0.5*(topology->cellLengths[1][loopA-1] + topology->cellLengths[1][loopA]);
    }
    fprintf(outFile,"%e\n",currYCoord);

    // Export Z Coordinates
    fprintf(outFile,"Z_COORDINATES %d double\n",(int)topology->cellLengths[2].size());
    currZCoord = topology->domainSizeMin[2];
    for(size_t loopA=1;loopA<topology->cellLengths[2].size();loopA++){
      fprintf(outFile,"%e ",currZCoord);
      currZCoord += 0.5*(topology->cellLengths[2][loopA-1] + topology->cellLengths[2][loopA]);
    }
    fprintf(outFile,"%e\n",currZCoord);
  }

  // Export Point quantities
  fprintf(outFile,"POINT_DATA %d\n",topology->totalCells);

  // Export Normal sign
  double* normSignX = new double[topology->totalCells];
  double* normSignY = new double[topology->totalCells];
  double* normSignZ = new double[topology->totalCells];
  int* counterVec = new int[topology->totalCells];
  for (int loopA=0;loopA<topology->totalCells;loopA++){
    normSignX[loopA] = 0.0;
    normSignY[loopA] = 0.0;
    normSignZ[loopA] = 0.0;
    counterVec[loopA] = 0;
  }

  int currCell = 0;
  for(size_t loopA=0;loopA<topology->faceConnections.size();loopA++){
    if(topology->faceCells[loopA].size() == 1){
      currCell = topology->faceCells[loopA][0];
      counterVec[currCell]++;
      normSignX[currCell] += topology->faceNormal[loopA][0];
      normSignY[currCell] += topology->faceNormal[loopA][1];
      normSignZ[currCell] += topology->faceNormal[loopA][2];
    }
  }

  for(int loopA=0;loopA<topology->totalCells;loopA++){
    if(counterVec[loopA] > 0){
      normSignX[loopA] /= (double)counterVec[loopA];
      normSignY[loopA] /= (double)counterVec[loopA];
      normSignZ[loopA] /= (double)counterVec[loopA];
    }
  }
  fprintf(outFile,"VECTORS faceNormals double\n");
  // Print velocity components
  for (int loopA=0;loopA<topology->totalCells;loopA++){
    fprintf(outFile,"%e %e %e\n",normSignX[loopA],normSignY[loopA],normSignZ[loopA]);
  }
  delete [] normSignX;
  delete [] normSignY;
  delete [] normSignZ;
  delete [] counterVec;

  // Print Scalar Concentration
  fprintf(outFile,"SCALARS concentration double\n");
  fprintf(outFile,"LOOKUP_TABLE default\n");
  // Print Concentrations
  for (int loopA=0;loopA<topology->totalCells;loopA++){
    fprintf(outFile,"%e\n",cells[loopA].concentration);
  }

  //Print velocity
  fprintf(outFile,"VECTORS velocity double\n");
  // Print velocity components
  for (int loopA=0;loopA<topology->totalCells;loopA++){
    fprintf(outFile,"%e %e %e\n",cells[loopA].velocity[0],cells[loopA].velocity[1],cells[loopA].velocity[2]);
    //fprintf(outFile,"%e %e %e\n",cellPoints[loopA].filteredVel[0],cellPoints[loopA].filteredVel[1],cellPoints[loopA].filteredVel[2]);
  }

  // Print Pressure Gradient
  if (hasPressureGradient){
    fprintf(outFile,"VECTORS PressureGrad float\n");
    // Print pressure Gradient
    for (int loopA=0;loopA<topology->totalCells;loopA++){
      fprintf(outFile,"%e %e %e\n",cells[loopA].pressGrad[0],cells[loopA].pressGrad[1],cells[loopA].pressGrad[2]);
    }
  }

  // Print Relative Pressure
  if (hasRelativePressure){
    fprintf(outFile,"SCALARS RelPressure double\n");
    fprintf(outFile,"LOOKUP_TABLE default\n");
    // Print Relative Pressure
    for (int loopA=0;loopA<topology->totalCells;loopA++){
      fprintf(outFile,"%e\n",cells[loopA].relPressure);
    }
  }

  // Print Relative Pressure
  if (hasRelativePressure){
    fprintf(outFile,"SCALARS GradientMonitor double\n");
    fprintf(outFile,"LOOKUP_TABLE default\n");
    // Print Relative Pressure
    for (int loopA=0;loopA<topology->totalCells;loopA++){
      fprintf(outFile,"%e\n",cells[loopA].auxVector[0]);
    }
  }

  // Print Reynolds Stresses
  if (hasReynoldsStress){
    fprintf(outFile,"TENSORS ReynoldsStress double\n");
    // Print Reynolds Stress Tensor
    for (int loopA=0;loopA<topology->totalCells;loopA++){
      fprintf(outFile,"%e %e %e\n",cells[loopA].ReStress[0],cells[loopA].ReStress[1],cells[loopA].ReStress[2]);
      fprintf(outFile,"%e %e %e\n",cells[loopA].ReStress[1],cells[loopA].ReStress[3],cells[loopA].ReStress[4]);
      fprintf(outFile,"%e %e %e\n",cells[loopA].ReStress[2],cells[loopA].ReStress[4],cells[loopA].ReStress[5]);
      fprintf(outFile,"\n");
    }
  }  

  // ==================
  // EXPORT DERIVATIVES
  // ==================

  // First and Second Derivatives
  MRIDoubleMat firstDerivs;
  MRIDoubleMat secondDerivs;
  firstDerivs.resize(kNumberOfDimensions);
  secondDerivs.resize(kNumberOfDimensions);
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    firstDerivs[loopA].resize(kNumberOfDimensions);
    secondDerivs[loopA].resize(kNumberOfDimensions);
  }

  // Print the Velocity Gradient
  fprintf(outFile,"TENSORS VelocityGradient double\n");
  // Print Reynolds Stress Tensor
  for (int loopA=0;loopA<topology->totalCells;loopA++){
    evalSpaceDerivs(loopA,threshold,firstDerivs,secondDerivs);
    fprintf(outFile,"%e %e %e\n",firstDerivs[0][0],firstDerivs[0][1],firstDerivs[0][2]);
    fprintf(outFile,"%e %e %e\n",firstDerivs[1][0],firstDerivs[1][1],firstDerivs[1][2]);
    fprintf(outFile,"%e %e %e\n",firstDerivs[2][0],firstDerivs[2][1],firstDerivs[2][2]);
    fprintf(outFile,"\n");
  }

  // Print the Velocity Gradient
  fprintf(outFile,"TENSORS VelocityCurvature double\n");
  // Print Reynolds Stress Tensor
  for (int loopA=0;loopA<topology->totalCells;loopA++){
    evalSpaceDerivs(loopA,threshold,firstDerivs,secondDerivs);
    fprintf(outFile,"%e %e %e\n",secondDerivs[0][0],secondDerivs[0][1],secondDerivs[0][2]);
    fprintf(outFile,"%e %e %e\n",secondDerivs[1][0],secondDerivs[1][1],secondDerivs[1][2]);
    fprintf(outFile,"%e %e %e\n",secondDerivs[2][0],secondDerivs[2][1],secondDerivs[2][2]);
    fprintf(outFile,"\n");
  }

  // EXPORT OUTPUTS
  for(size_t loopA=0;loopA<outputs.size();loopA++){
    // Print Header
    if(outputs[loopA].totComponents == 1){
      fprintf(outFile,"%s\n",string("SCALARS " + outputs[loopA].name + " double").c_str());
      fprintf(outFile,"LOOKUP_TABLE default\n");
      for (int loopB=0;loopB<topology->totalCells;loopB++){
        fprintf(outFile,"%e\n",outputs[loopA].values[loopB]);
      }
    }else{
      int count = 0;
      fprintf(outFile,"%s\n",string("VECTORS " + outputs[loopA].name + " double").c_str());
      for (int loopB=0;loopB<topology->totalCells;loopB++){
        for(int loopC=0;loopC<outputs[loopA].totComponents;loopC++){
          fprintf(outFile,"%e ",outputs[loopA].values[count]);
          count++;
        }
        fprintf(outFile,"\n");
      }
    }
  }

  // Print Tagging
  if(cellTags.size() > 0){
    fprintf(outFile,"SCALARS Tags double\n");
    fprintf(outFile,"LOOKUP_TABLE default\n");
    for (int loopA=0;loopA<topology->totalCells;loopA++){
      fprintf(outFile,"%d\n",cellTags[loopA]);
    }
  }

  // Close File
  fclose(outFile);
}

// Eval total Number of Vortices
int MRIScan::evalTotalVortex(){
  // Init
  int total = 0;
  int totalSlices = 0;
  int totalStars = 0;
  // YZ Planes
  totalSlices = topology->cellTotals[0];
  totalStars = (topology->cellTotals[1]+1)*(topology->cellTotals[2]+1);
  total += totalSlices * totalStars;
  // XZ Planes
  totalSlices = topology->cellTotals[1];
  totalStars = (topology->cellTotals[0]+1)*(topology->cellTotals[2]+1);
  total += totalSlices * totalStars;
  // XY Planes
  totalSlices = topology->cellTotals[2];
  totalStars = (topology->cellTotals[0]+1)*(topology->cellTotals[1]+1);
  total += totalSlices * totalStars;
  // Return Value
  return total;
}

// ======================
// Read List of Row Files
// ======================
void MRIScan::readRAWFileSequence(std::string fileListName){

  // Init
  MRIImageData data;
  std::vector<std::string> fileList;

  // Read File List
  MRIUtils::readFileList(fileListName,fileList);

  // Read the first File
  std::string currFile = fileList[0];

  // Read
  readRawImage(currFile,data);

  // Set Totals
  topology->cellTotals[0] = data.sizeX;
  topology->cellTotals[1] = data.sizeY;
  topology->cellTotals[2] = fileList.size();

  // Set cell Lengths 1.0 Uniform in all directions
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    for(size_t loopB=0;loopB<topology->cellLengths[loopA].size();loopB++){
      topology->cellLengths[loopA][loopB] = 1.0;
    }
  }

  // Set total Points
  topology->totalCells = topology->cellTotals[0] * topology->cellTotals[1] * topology->cellTotals[2];

  // Intialize Cells
  MRICell myCellPoint;
  myCellPoint.concentration = 0.0;
  myCellPoint.velocity[0] = 0.0;
  myCellPoint.velocity[1] = 0.0;
  myCellPoint.velocity[2] = 0.0;
  // Resize: CHECK!!!
  cells.resize(topology->totalCells,myCellPoint);

  // Fill Concentration with First Image Data
  for(int loopA=0;loopA<data.sizeX*data.sizeY;loopA++){
    cells[loopA].concentration = data.rawData[loopA];
  }

  // Intialize how many cells read so far
  int readSoFar = data.sizeX*data.sizeY;
  // Loop through all other images
  for(unsigned int loopA=1;loopA<fileList.size();loopA++){
    currFile = fileList[loopA];
    // Read
    readRawImage(currFile,data);
    // Check Consistency
    if((data.sizeX != topology->cellTotals[0])||(data.sizeY != topology->cellTotals[1])){
      throw new MRIException("Error: Image Sequence not Consistent!");
    }
    // Fill Concentration with First Image Data
    for(int loopB=readSoFar;loopB<readSoFar + data.sizeX*data.sizeY;loopB++){
      cells[loopB].concentration = data.rawData[loopB];
    }
    // Increment readSoFar
    readSoFar += data.sizeX*data.sizeY;
  }
}

// ====================
// Read Raw File Header
// ====================
void readRawFileHeader(int &sizeX,int &sizeY,int &numberOfBytes,FILE *fp){
  // Read Dimensions
  int wordCount = 0;
  int currentByte = 0;
  size_t fres;
  std::string fileWord = "";
  while(wordCount<4){
    currentByte = 0;
    fileWord = "";
    while((currentByte != 10)&&(currentByte != 32)){
      fres = fread(&currentByte, 1, 1, fp);
      if((currentByte != 10)&&(currentByte != 32)){
        fileWord = fileWord + char(currentByte);
      }
    }
    // Recover Current Word
    if(wordCount == 1){
      // Size X
      sizeX = atoi(fileWord.c_str());
    }else if (wordCount == 2){
      // Size Y
      sizeY = atoi(fileWord.c_str());
    }else if (wordCount == 3){
      // Number Of Bytes: CHECK!!!
      //numberOfBytes = atoi(fileWord.c_str());
      //numberOfBytes = (log(numberOfBytes+1)/log(2.0));
    }

    // Update Value
    wordCount++;
  }
}


// ===================
// Read Raw Image Data
// ===================
int MRIScan::readRawImage(std::string FileName, MRIImageData &data){
  // Open File
  FILE *fp = NULL;
  fp = fopen(FileName.c_str(),"rb");
  if (fp == NULL){
    printf("Error: Failed to read volume  %s\n",FileName.c_str());
    return -1;
  }
  // Read Raw File Header
  int numberOfBytes = 0;
  readRawFileHeader(data.sizeX,data.sizeY,numberOfBytes,fp);

  // Set the size of the voxel information
  int size = data.sizeX * data.sizeY;

  // Read voxel information
  data.rawData = new short[size];
  int len = fread(data.rawData,1,size,fp);
  if (len != size){
    printf("Error: Failed to read enough data !\n");
  }

  // Close File
  fclose(fp);

  // Return
  return 0;
}

// ===================
// READ EXPANSION FILE
// ===================
void readExpansionFile(std::string fileName,int* tot,
                       std::vector<double> &lengthX,
                       std::vector<double> &lengthY,
                       std::vector<double> &lengthZ,
                       double* minlimits,double* maxlimits,MRIExpansion* &exp){

  // ASSIGN FILE
  int lineCount = 0;
  std::vector<std::string> ResultArray;
  std::string Buffer;
  std::ifstream inFile;
  inFile.open(fileName.c_str());

  // GET TOTAL CELLS
  lineCount++;
  std::getline(inFile,Buffer);
  ResultArray = MRIUtils::extractSubStringFromBufferMS(Buffer);
  tot[0] = atoi(ResultArray[0].c_str());
  tot[1] = atoi(ResultArray[1].c_str());
  tot[2] = atoi(ResultArray[2].c_str());

  printf("TOTALS %d %d %d\n",tot[0],tot[1],tot[2]);

  // GET CELL X LENGTHS
  int lengthCount = 0;
  while(lengthCount<tot[0]){
    std::getline(inFile,Buffer);
    boost::trim(Buffer);
    ResultArray = MRIUtils::extractSubStringFromBufferMS(Buffer);
    for(size_t loopA=0;loopA<ResultArray.size();loopA++){
      // Assign Length
      lengthX.push_back(atof(ResultArray[loopA].c_str()));
      // Update Counter
      lengthCount++;
    }
  }

  // GET CELL Y LENGTHS
  lengthCount = 0;
  while(lengthCount<tot[1]){
    std::getline(inFile,Buffer);
    boost::trim(Buffer);
    ResultArray = MRIUtils::extractSubStringFromBufferMS(Buffer);
    for(size_t loopA=0;loopA<ResultArray.size();loopA++){
      // Assign Length
      lengthY.push_back(atof(ResultArray[loopA].c_str()));
      // Update Counter
      lengthCount++;
    }
  }

  // GET CELL Z LENGTHS
  lengthCount = 0;
  while(lengthCount<tot[2]){
    std::getline(inFile,Buffer);
    boost::trim(Buffer);
    ResultArray = MRIUtils::extractSubStringFromBufferMS(Buffer);
    for(size_t loopA=0;loopA<ResultArray.size();loopA++){
      // Assign Length
      lengthZ.push_back(atof(ResultArray[loopA].c_str()));
      // Update Counter
      lengthCount++;
    }
  }

  // GET MIN LIMITS
  lineCount++;
  std::getline(inFile,Buffer);
  ResultArray = MRIUtils::extractSubStringFromBufferMS(Buffer);
  minlimits[0] = atof(ResultArray[0].c_str());
  minlimits[1] = atof(ResultArray[1].c_str());
  minlimits[2] = atof(ResultArray[2].c_str());

  // GET MAX LIMITS
  lineCount++;
  std::getline(inFile,Buffer);
  ResultArray = MRIUtils::extractSubStringFromBufferMS(Buffer);
  maxlimits[0] = atof(ResultArray[0].c_str());
  maxlimits[1] = atof(ResultArray[1].c_str());
  maxlimits[2] = atof(ResultArray[2].c_str());

  // GET EXPANSION COEFFICIENTS
  std::vector<double> tempExpansion;
  while(std::getline(inFile,Buffer)){
    ResultArray = MRIUtils::extractSubStringFromBufferMS(Buffer);
    tempExpansion.push_back(atof(ResultArray[0].c_str()));
  }

  // CREATE EXPANSION
  exp = new MRIExpansion((int)tempExpansion.size()-3);
  exp->FillFromVector(tempExpansion);

  // CLOSE FILE
  inFile.close();
}

// =============================
// READ SCAN FROM EXPANSION FILE
// =============================
void MRIScan::readFromExpansionFile(std::string fileName,bool applyThreshold, int thresholdType,double thresholdRatio){

  // Allocate Variables
  int tot[3];
  std::vector<double> lengthX;
  std::vector<double> lengthY;
  std::vector<double> lengthZ;
  double minlimits[3];
  double maxlimits[3];
  MRIExpansion* exp = NULL;

  // Read Quantities From File
  readExpansionFile(fileName,tot,lengthX,lengthY,lengthZ,minlimits,maxlimits,exp);

  // SET UP SCAN QUANTITIES
  // CELL TOTALS
  topology->cellTotals[0] = tot[0];
  topology->cellTotals[1] = tot[1];
  topology->cellTotals[2] = tot[2];

  // CELL LENGTHS
  topology->cellLengths.resize(kNumberOfDimensions);
  // X
  for(size_t loopA=0;loopA<lengthX.size();loopA++){
    topology->cellLengths[0].push_back(lengthX[loopA]);
  }
  // Y
  for(size_t loopA=0;loopA<lengthY.size();loopA++){
    topology->cellLengths[1].push_back(lengthY[loopA]);
  }
  // Z
  for(size_t loopA=0;loopA<lengthZ.size();loopA++){
    topology->cellLengths[2].push_back(lengthZ[loopA]);
  }

  // DIMENSIONS
  // MIN
  topology->domainSizeMin[0] = minlimits[0];
  topology->domainSizeMin[1] = minlimits[1];
  topology->domainSizeMin[2] = minlimits[2];
  // MAX
  topology->domainSizeMax[0] = maxlimits[0];
  topology->domainSizeMax[1] = maxlimits[1];
  topology->domainSizeMax[2] = maxlimits[2];
  // MRI EXPANSION
  expansion = new MRIExpansion(exp);

  // APPLY THRESHOLD TO EXPANSION
  if(applyThreshold){
    expansion->ApplyVortexThreshold(thresholdType,thresholdRatio);
  }

  // INITIALIZE VALUES
  hasPressureGradient = false;
  hasRelativePressure = false;
  scanTime = 0.0;
  maxVelModule = 0.0;

  // INITIALIZE SCAN
  topology->totalCells = topology->cellTotals[0]*topology->cellTotals[1]*topology->cellTotals[2];
  // CREATE NEW CELL
  MRICell newCell;
  // INITIALIZE QUANTITIES TO ZERO
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // ADD IT TO CELL POINTS
    cells.push_back(newCell);
  }

  // INITIALIZE POSITIONS
  int intCoords[3] = {0};
  double Pos[3] = {0.0};
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    mapIndexToCoords(loopA,intCoords);
    mapCoordsToPosition(intCoords,true,Pos);
    cells[loopA].setQuantity(kQtyPositionX,topology->domainSizeMin[0] + Pos[0]);
    cells[loopA].setQuantity(kQtyPositionY,topology->domainSizeMin[1] + Pos[1]);
    cells[loopA].setQuantity(kQtyPositionZ,topology->domainSizeMin[2] + Pos[2]);
  }

  // REORDER MODEL
  reorderScan();

  // REBUILD SCAN
  //RebuildFromExpansion(expansion,true);
  // No Constant Flux
  rebuildFromExpansion(expansion,false);
}

// ====================
// WRITE EXPANSION FILE
// ====================
void MRIScan::writeExpansionFile(std::string fileName){
  // Open Output File
  FILE* fid;
  fid = fopen(fileName.c_str(),"w");

  // WRITE TOTAL CELLS
  fprintf(fid,"%15d %15d %15d\n",topology->cellTotals[0],topology->cellTotals[1],topology->cellTotals[2]);


  // WRITE CELL LENGTHS X
  for(int loopA=0;loopA<topology->cellTotals[0];loopA++){
    fprintf(fid,"%15.6e\n",topology->cellLengths[0][loopA]);
  }

  // WRITE CELL LENGTHS Y
  for(int loopA=0;loopA<topology->cellTotals[1];loopA++){
    fprintf(fid,"%15.6e\n",topology->cellLengths[1][loopA]);
  }

  // WRITE CELL LENGTHS Z
  for(int loopA=0;loopA<topology->cellTotals[2];loopA++){
    fprintf(fid,"%15.6e\n",topology->cellLengths[2][loopA]);
  }

  // MIN DOMAIN SIZE
  fprintf(fid,"%15.6e %15.6e %15.6e\n",topology->domainSizeMin[0],topology->domainSizeMin[1],topology->domainSizeMin[2]);
  // MAX DOMAIN SIZE
  fprintf(fid,"%15.6e %15.6e %15.6e\n",topology->domainSizeMax[0],topology->domainSizeMax[1],topology->domainSizeMax[2]);

  // WRITE EXPANSION COEFFICIENTS
  // Write Constant Flux Components
  for(int loopA=0;loopA<3;loopA++){
    fprintf(fid,"%15.6e\n",expansion->constantFluxCoeff[loopA]);
  }
  // Write Vortex Component
  for(int loopA=0;loopA<expansion->totalVortices;loopA++){
    fprintf(fid,"%15.6e\n",expansion->vortexCoeff[loopA]);
  }
  // Close Output file
  fclose(fid);
}

// ==================
// THRESHOLD QUANTITY
// ==================
void MRIScan::thresholdQuantity(int qtyID,double threshold){
  // DECLARE
  WriteSchMessage(std::string("Applying threshold...\n"));
  double centerCellValue = 0.0;
  // LOOP ON ALL CELLS
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Get Value in Current Cell
    centerCellValue = cells[loopA].getQuantity(qtyID);
    // Apply Threshold
    if (centerCellValue<threshold){
      cells[loopA].setQuantity(qtyID,0.0);
    }
  }
}

// ====================================================================
// DETERMINE THREE-DIMENSIONAL COMPONENTS OF THE EXPANSION COEFFICIENTS
// ====================================================================
void MRIScan::evalSMPVortexCriteria(MRIExpansion* exp){
  // LOOP ON CELLS
  MRIIntVec idx;
  MRIOutput out1("SMPVortexCriterion",3);
  double avVortexIndex = 0.0;
  double totalIntensity = 0.0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Loop on the dimensions
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      // Determine the idexes for the adjacent vortices
      getNeighborVortexes(loopA,loopB,idx);
      // Average the values
      avVortexIndex = 0.25*(exp->vortexCoeff[idx[0]] +
                            exp->vortexCoeff[idx[1]] +
                            exp->vortexCoeff[idx[2]] +
                            exp->vortexCoeff[idx[3]]);
      // Get total Intensity
      totalIntensity += (exp->vortexCoeff[idx[0]]*exp->vortexCoeff[idx[0]] +
                         exp->vortexCoeff[idx[1]]*exp->vortexCoeff[idx[1]] +
                         exp->vortexCoeff[idx[2]]*exp->vortexCoeff[idx[2]] +
                         exp->vortexCoeff[idx[3]]*exp->vortexCoeff[idx[3]]);
      // Assign To cell
      out1.values.push_back(avVortexIndex);
    }
  }
  // Add output quantities
  outputs.push_back(out1);
  // Printout
  printf("Total SMP Intensities: %e\n",totalIntensity);
}

// ===============================================
// GET TOTAL NUMBER OF FACES FOR STRUCTURED LAYOUT
// ===============================================
int MRIScan::getTotalFaces(){
  return topology->cellTotals[0]*topology->cellTotals[1]*(topology->cellTotals[2] + 1) +
         topology->cellTotals[1]*topology->cellTotals[2]*(topology->cellTotals[0] + 1) +
         topology->cellTotals[2]*topology->cellTotals[0]*(topology->cellTotals[1] + 1);
}

// Get Face Starting From Cell and Unit Vector
int MRIScan::getFacewithCellVector(int CurrentCell, double* UnitVector){
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
        throw new MRIException("Internal Error: Problems in GetFacewithCellVector");
    }
  // Get Adj Face
  resultFace = getAdjacentFace(CurrentCell,AdjType);
  return resultFace;
}

// ==========================
// GET LOCAL FACE CONNECTIONS
// ==========================
void getFaceConnections(int faceID, std::vector<int> cellConnections, std::vector<int> &faceIds){
  faceIds.clear();
  switch(faceID){
    case 0:
      faceIds.push_back(cellConnections[0]);
      faceIds.push_back(cellConnections[1]);
      faceIds.push_back(cellConnections[3]);
      faceIds.push_back(cellConnections[2]);
      break;
    case 1:
      faceIds.push_back(cellConnections[4]);
      faceIds.push_back(cellConnections[6]);
      faceIds.push_back(cellConnections[7]);
      faceIds.push_back(cellConnections[5]);
      break;
    case 2:
      faceIds.push_back(cellConnections[0]);
      faceIds.push_back(cellConnections[2]);
      faceIds.push_back(cellConnections[6]);
      faceIds.push_back(cellConnections[4]);
      break;
    case 3:
      faceIds.push_back(cellConnections[1]);
      faceIds.push_back(cellConnections[5]);
      faceIds.push_back(cellConnections[7]);
      faceIds.push_back(cellConnections[3]);
      break;
    case 4:
      faceIds.push_back(cellConnections[0]);
      faceIds.push_back(cellConnections[4]);
      faceIds.push_back(cellConnections[5]);
      faceIds.push_back(cellConnections[1]);
      break;
    case 5:
      faceIds.push_back(cellConnections[2]);
      faceIds.push_back(cellConnections[3]);
      faceIds.push_back(cellConnections[7]);
      faceIds.push_back(cellConnections[6]);
      break;
  }
}

// ====================================
// GET LOCAL EDGE CONNECTIONS FROM FACE
// ====================================
void getEdgeConnections(int EdgeID, std::vector<int> faceConnections, int* edgeIds){
    switch(EdgeID){
      case 0:
        edgeIds[0] = faceConnections[0];
        edgeIds[1] = faceConnections[1];
        break;
      case 1:
        edgeIds[0] = faceConnections[1];
        edgeIds[1] = faceConnections[2];
        break;
      case 2:
        edgeIds[0] = faceConnections[2];
        edgeIds[1] = faceConnections[3];
        break;
      case 3:
        edgeIds[0] = faceConnections[3];
        edgeIds[1] = faceConnections[0];
        break;
    }
}

// ======================
// GET VORTEX COEFFICIENT
// ======================
double MRIScan::getEdgeFaceVortexCoeff(int edgeID, int faceID){
  double edgeDirVector[3];
  double edgeFaceVector[3];
  double resVec[3];
  double res = 0.0;
  // Get Vectors
  getEdgeDirection(edgeID,edgeDirVector);
  getEdgeToFaceDirection(edgeID,faceID,edgeFaceVector);

  // Eval Vector Product
  resVec[0] = edgeDirVector[1] * edgeFaceVector[2] - edgeFaceVector[1] * edgeDirVector[2];
  resVec[1] = edgeDirVector[2] * edgeFaceVector[0] - edgeFaceVector[2] * edgeDirVector[0];
  resVec[2] = edgeDirVector[0] * edgeFaceVector[1] - edgeFaceVector[0] * edgeDirVector[1];
  double modulus = (resVec[0] * resVec[0] + resVec[1] * resVec[1] + resVec[2] * resVec[2]);
  resVec[0] = resVec[0]/modulus;
  resVec[1] = resVec[1]/modulus;
  resVec[2] = resVec[2]/modulus;
  // Get Sign
  res = resVec[0] * topology->faceNormal[faceID][0] + 
        resVec[1] * topology->faceNormal[faceID][1] + 
        resVec[2] * topology->faceNormal[faceID][2];
  return round(res);
}

// ==========================
// GET EDGE TO FACE DIRECTION
// ==========================
void MRIScan::getEdgeToFaceDirection(int edgeID, int faceID, double* edgeFaceVector){
  // Declare
  double ec[3] = {0.0};
  double fc[3] = {0.0};

  // Get Edge Center
  getEdgeCenter(edgeID,ec);

  // Get Face Center
  getFaceCenter(faceID,fc);

  // Get the versor
  edgeFaceVector[0] = fc[0] - ec[0];
  edgeFaceVector[1] = fc[1] - ec[1];
  edgeFaceVector[2] = fc[2] - ec[2];
  double modulus = (edgeFaceVector[0] * edgeFaceVector[0] + edgeFaceVector[1] * edgeFaceVector[1] + edgeFaceVector[2] * edgeFaceVector[2]);
  edgeFaceVector[0] = edgeFaceVector[0]/modulus;
  edgeFaceVector[1] = edgeFaceVector[1]/modulus;
  edgeFaceVector[2] = edgeFaceVector[2]/modulus;
}

// ===============
// GET EDGE CENTER
// ===============
void MRIScan::getEdgeCenter(int edgeID, double* ec){
  int node1 = 0;
  int node2 = 0;
  double node1Pos[3] = {0.0};
  double node2Pos[3] = {0.0};

  // Get The Two Nodes
  node1 = topology->edgeConnections[edgeID][0];
  node2 = topology->edgeConnections[edgeID][1];

  // Eval Auxiliary Node Coordinates
  node1Pos[0] = topology->auxNodesCoords[node1][0];
  node1Pos[1] = topology->auxNodesCoords[node1][1];
  node1Pos[2] = topology->auxNodesCoords[node1][2];
  node2Pos[0] = topology->auxNodesCoords[node2][0];
  node2Pos[1] = topology->auxNodesCoords[node2][1];
  node2Pos[2] = topology->auxNodesCoords[node2][2];

  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    ec[loopA] = 0.5*(node1Pos[loopA] + node2Pos[loopA]);
  }
}

// ===============
// GET FACE CENTER
// ===============
void MRIScan::getFaceCenter(int faceID, double* fc){
  int currNode = 0;
  double pos[3] = {0.0};

  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    fc[loopA] = 0.0;
  }

  for(size_t loopA=0;loopA<topology->faceConnections[faceID].size();loopA++){
    currNode = topology->faceConnections[faceID][loopA];
    pos[0] = topology->auxNodesCoords[currNode][0];
    pos[1] = topology->auxNodesCoords[currNode][1];
    pos[2] = topology->auxNodesCoords[currNode][2];
    for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
      fc[loopA] += pos[loopA];
    }
  }
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    fc[loopA] /= (double)topology->faceConnections[faceID].size();
  }
}

// ========================================================
// EVAL THE AVERGAGE ERROR BETWEEN FILTEREDVEL AND VELOCITY
// ========================================================
void MRIScan::recoverGlobalErrorEstimates(double& AvNormError, double& AvAngleError){
  // Init
  AvNormError = 0.0;
  AvAngleError = 0.0;
  double diffNorm;
  double cosAlpha,currentAlpha;
  MRIDoubleVec diffVel(3,0.0);
  MRIDoubleVec normVel(3,0.0);
  MRIDoubleVec normFilterVel(3,0.0);
  // Loop
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Eval Velocity Difference
    diffVel[0] = cells[loopA].velocity[0]-cells[loopA].auxVector[0];
    diffVel[1] = cells[loopA].velocity[1]-cells[loopA].auxVector[1];
    diffVel[2] = cells[loopA].velocity[2]-cells[loopA].auxVector[2];
    diffNorm = MRIUtils::do3DEucNorm(diffVel);
    AvNormError = AvNormError + diffNorm;
    // Eval Velocity Angle
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      normVel[loopB] = cells[loopA].velocity[loopB];
      normFilterVel[loopB] = cells[loopA].auxVector[loopB];
    }
    // Normalize
    MRIUtils::normalize3DVector(normVel);
    MRIUtils::normalize3DVector(normFilterVel);
    cosAlpha = normVel[0]*normFilterVel[0]+
               normVel[1]*normFilterVel[1]+
               normVel[2]*normFilterVel[2];
    if(cosAlpha>1.0) cosAlpha = 1.0;
    if(cosAlpha<-1.0) cosAlpha = -1.0;
    currentAlpha = acos(cosAlpha);
    AvAngleError = AvAngleError + currentAlpha;
  }
  // Eval Average Values
  if (topology->totalCells>0){
    AvNormError = (AvNormError/topology->totalCells);
    AvAngleError = (AvAngleError/topology->totalCells);
  }else{
    AvNormError = -1.0;
    AvAngleError = -1.0;
  }
}



// ===================================
// BUILD TOPOLOGY VECTORS FOR PARMETIS
// ===================================
void MRIScan::buildMetisConnectivities(int *eptr,int *eind){
  int currNode = 0;
  // Allocate the pointer
  eptr = new int(topology->cellConnections.size());

  // Initialize Counter
  int count = 0;
  // Loop through the cells
  for(size_t loopA=0;loopA<topology->cellConnections.size();loopA++){
    // Assign the pointer
    eptr[loopA] = count;
    for(size_t loopB=0;loopB<topology->cellConnections[loopA].size();loopB++){
      // Increment Counter
      count++;
    }
  }
  // Initialize the Connectivities
  eind = new int(count);

  // Loop through the cells
  count = 0;
  for(size_t loopA=0;loopA<topology->cellConnections.size();loopA++){
    for(size_t loopB=0;loopB<topology->cellConnections[loopA].size();loopB++){
      // Get Current Node Number
      currNode = topology->cellConnections[loopA][loopB];
      // Increment Counter
      count++;
      // Assign Node
      eind[count] = currNode;
    }
  }
}

// Print VTK OPTIONS
void printVTKOptions(vtkStructuredPointsOptionRecord opts){
  printf("\n");
  if(opts.isASCII){
    printf("ASCII FILE FORMAT\n");
  }else{
    printf("BINARY FILE FORMAT\n");
  }
  if(opts.isValidDataset){
    printf("STRUCTURED POINTS\n");
  }else{
    printf("OTHER VTK FORMAT\n");
  }
  printf("GRID SIZE: %d %d %d\n",opts.dimensions[0],opts.dimensions[1],opts.dimensions[2]);
  printf("GRID ORIGIN: %e %e %e\n",opts.origin[0],opts.origin[1],opts.origin[2]);
  printf("GRID SPACING: %e %e %e\n",opts.spacing[0],opts.spacing[1],opts.spacing[2]);
  printf("\n");
}

// INIT OPTIONS
void initVTKStructuredPointsOptions(vtkStructuredPointsOptionRecord &opts){
  opts.isASCII = false;
  opts.isValidDataset = false;
  opts.numDefined = 5;
  // Size
  opts.dimensions[0] = 0;
  opts.dimensions[1] = 0;
  opts.dimensions[2] = 0;
  // Origin
  opts.origin[0] = 0.0;
  opts.origin[1] = 0.0;
  opts.origin[2] = 0.0;
  // Spacing
  opts.spacing[0] = 0.0;
  opts.spacing[1] = 0.0;
  opts.spacing[2] = 0.0;
  // Is Defined
  for(int loopA=0;loopA<opts.numDefined;loopA++){
    opts.isDefined[loopA] = false;
  }
}

// Read Scalars From VTK Legacy
void readVTKScalar(ifstream& vtkFile, int& lineCount,int linesToRead,MRIDoubleVec& vtkScalar){
  vtkScalar.clear();
  string buffer;
  double value;
  vector<string> tokenizedString;
  // Skip first two lines
  std::getline(vtkFile,buffer);
  //printf("Skipping; %s\n",buffer.c_str());
  lineCount++;
  std::getline(vtkFile,buffer);
  //printf("Skipping; %s\n",buffer.c_str());
  lineCount++;
  for(int loopA=0;loopA<linesToRead;loopA++){
    std::getline(vtkFile,buffer);
    //printf("Reading: %s\n",buffer.c_str());
    boost::trim(buffer);
    lineCount++;
    if(!buffer.empty()){
      boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
      for(size_t loopB=0;loopB<tokenizedString.size();loopB++){
        try{
          value = atof(tokenizedString[loopB].c_str());
          vtkScalar.push_back(value);
        }catch(...){
          return;
        }
      }
    }else{
      break;
    }
  }
}

// Read Vectors From VTK Legacy
void readVTKVector(ifstream& vtkFile, int& lineCount,int linesToRead,MRIDoubleMat& vtkVector){
  vtkVector.clear();
  MRIDoubleVec temp;
  MRIDoubleVec store;
  string buffer;
  double value;
  vector<string> tokenizedString;
  // Skip first line
  std::getline(vtkFile,buffer);
  lineCount++;
  for(int loopA=0;loopA<linesToRead;loopA++){
    std::getline(vtkFile,buffer);
    boost::trim(buffer);
    lineCount++;
    if(!buffer.empty()){
      boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
      for(size_t loopB=0;loopB<tokenizedString.size();loopB++){
        try{
          value = atof(tokenizedString[loopB].c_str());
          temp.push_back(value);
        }catch(...){
          //printf("Exiting with: %s\n",buffer.c_str());
          break;
        }
      }
    }else{
      //printf("Exiting with: %s\n",buffer.c_str());
      break;
    }
  }
  // Fill Matrix
  for(int loopA=0;loopA<temp.size()/3;loopA++){
    store.clear();
    for(int loopB=0;loopB<3;loopB++){
    }
    vtkVector.push_back(store);
  }
}

// ==========================
// READ VTK STRUCTURED POINTS
// ==========================
void MRIScan::readVTKStructuredPoints(std::string vtkFileName, bool DoReorderCells){

  // Init totalCellPoints
  topology->totalCells = 0;

  // Assign File
  WriteSchMessage(std::string("\n"));
  WriteSchMessage(std::string("--- READING STRUCTURED POINT FILE\n"));
  WriteSchMessage(std::string("\n"));
  std::ifstream vtkFile;
  WriteSchMessage(std::string("Open File: ") + vtkFileName + std::string("\n"));
  vtkFile.open(vtkFileName.c_str());

  // Create and initialize vtkOption Vector
  vtkStructuredPointsOptionRecord vtkOptions;
  initVTKStructuredPointsOptions(vtkOptions);

  // Read Through and look for options
  std::vector<std::string> tokenizedString;
  WriteSchMessage(std::string("Computing input file size..."));
  int totalLinesInFile = 0;
  std::string Buffer;
  while (std::getline(vtkFile,Buffer)){
    boost::split(tokenizedString, Buffer, boost::is_any_of(" ,"), boost::token_compress_on);
    // Check if you find options
    assignVTKOptions(totalLinesInFile,tokenizedString, vtkOptions);
    // Increase line number
    totalLinesInFile++;
  }
  vtkOptions.dataBlockStart.push_back(totalLinesInFile);
  vtkOptions.dataBlockType.push_back(0);
  vtkOptions.dataBlockRead.push_back(false);

  // Done: Computing Input File Size
  WriteSchMessage(std::string("Done.\n"));

  // Check if all properties were defined
  bool fileOK = true;
  for(int loopA=0;loopA<vtkOptions.numDefined;loopA++){
    fileOK = fileOK && vtkOptions.isDefined[loopA];
    if(!fileOK){
      printf("ERROR: DEFINITION %d\n",loopA);
    }
  }
  if(!fileOK){
    WriteSchMessage(std::string("ERROR: Invalid VTK File format.\n"));
    WriteSchMessage(std::string("\n"));
    exit(1);
  }

  // Print VTK OPTIONS
  printVTKOptions(vtkOptions);

  // Creating Grid Geometry from Options
  WriteSchMessage(std::string("Creating Grid Geometry ..."));
  createGridFromVTKStructuredPoints(vtkOptions);
  WriteSchMessage(std::string("Done.\n"));

  // Reset File
  vtkFile.clear();
  vtkFile.seekg(0, std::ios::beg);

  // Read Concentration and Velocities
  MRIDoubleVec vtkScalar;
  MRIDoubleMat vtkVector;
  int linesToRead = 0;
  totalLinesInFile = 0;
  for(int loopA=0;loopA<vtkOptions.dataBlockStart.size();loopA++){
    if(vtkOptions.dataBlockRead[loopA]){
      // Go to next block start
      while(totalLinesInFile<vtkOptions.dataBlockStart[loopA]){
        std::getline(vtkFile,Buffer);
        totalLinesInFile++;
      }
      // Read Scalars Or Vectors
      if(vtkOptions.dataBlockType[loopA] == 0){
        // Read Scalars
        WriteSchMessage(std::string("Reading Scalars...\n"));
        vtkScalar.clear();
        linesToRead = vtkOptions.dataBlockStart[loopA + 1] - vtkOptions.dataBlockStart[loopA] - 2;
        readVTKScalar(vtkFile,totalLinesInFile,linesToRead,vtkScalar);
        if(vtkScalar.size() != topology->totalCells){
          printf("Scalar Size %d, Total Cells %d\n",(int)vtkScalar.size(), topology->totalCells);
          throw MRIException("ERROR: Total number of scalars differ from number of cells.\n");
        }
      }else{
        // Read Vectors
        WriteSchMessage(std::string("Reading Vectors...\n"));
        vtkVector.clear();
        linesToRead = vtkOptions.dataBlockStart[loopA + 1] - vtkOptions.dataBlockStart[loopA];
        readVTKVector(vtkFile,totalLinesInFile,linesToRead,vtkVector);
        if(vtkVector.size() != topology->totalCells){
          printf("Vector Size %d, Total Cells %d\n",(int)vtkVector.size(), topology->totalCells);
          throw MRIException("ERROR: Total number of vectors differs from number of cells.\n");
        }
      }
    }
  }

  // Transfer Scalars and Vectors to Cells
  // Scalars
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    cells[loopA].concentration = vtkScalar[loopA];
  }
  // Vectors
  maxVelModule = 0.0;
  double currModulus = 0.0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
      cells[loopA].velocity[0] = vtkVector[loopA][0];
      cells[loopA].velocity[1] = vtkVector[loopA][1];
      cells[loopA].velocity[2] = vtkVector[loopA][2];
      currModulus = sqrt((cells[loopA].velocity[0] * cells[loopA].velocity[0]) +
                         (cells[loopA].velocity[1] * cells[loopA].velocity[1]) +
                         (cells[loopA].velocity[2] * cells[loopA].velocity[2]));
      if(currModulus > maxVelModule){
        maxVelModule = currModulus;
      }
  }

  // Finished Reading File
  WriteSchMessage(std::string("File reading completed.\n"));

  // Close File
  vtkFile.close();

  // REORDER CELLS
  //if (DoReorderCells){
  //  ReorderScan();
  //}

  // WRITE STATISTICS
  std::string CurrentStats = WriteStatistics();
  WriteSchMessage(CurrentStats);

}

// ============================
// GET FACE ID FROM CELL NUMBER
// ============================
int MRIScan::getCellFaceID(int CellId,int FaceId){
  int count = 0;
  bool found = false;
  while((!found)&&(count<topology->cellFaces[CellId].size())){
    // Check if found
    found = (topology->cellFaces[CellId][count] == FaceId);
    // Update
    if(!found){
      count++;
    }
  }
  if(!found){
    throw MRIException("ERROR: Face not found in MRIScan::GetCellFaceID.\n");
  }
  // Return value
  return count;
}

// =======================================
// SET TO ZERO THE FACES NOT ON THE BORDER
// =======================================
void MRIScan::setWallFluxesToZero(bool* isFaceOnWalls, MRIDoubleVec& poissonSourceFaceVec){
  for(int loopA=0;loopA<topology->faceConnections.size();loopA++){
    if((topology->faceCells[loopA].size() > 1)&&(isFaceOnWalls[loopA])){
      poissonSourceFaceVec[loopA] = 0.0;
    }
  }
}

// ===========================================================
// COMPUTE TURBULENT VISCOSITY USING PANDTL MIXED LENGTH MODEL
// ===========================================================
void MRIScan::evalPradtlTurbViscosity(MRIDoubleMat cellDistance, MRIThresholdCriteria* threshold, double density, MRIDoubleMat& turbViscosity){
  int cellCount = 0;
  double sTerm  = 0.0;
  double qty    = 0.0;
  double currDist = 0.0;
  MRIDoubleVec tmp;
  
  // First and Second Derivatives
  MRIDoubleMat firstDerivs;
  MRIDoubleMat secondDerivs;
  firstDerivs.resize(kNumberOfDimensions);
  secondDerivs.resize(kNumberOfDimensions);
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    firstDerivs[loopA].resize(kNumberOfDimensions);
    secondDerivs[loopA].resize(kNumberOfDimensions);
  }
  
  turbViscosity.clear();
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Only Cells with significant Concentration
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){
      // Get Current Distance
      currDist = cellDistance[cellCount][0];
      // Eva Spatial Derivatives
      evalSpaceDerivs(loopA, threshold, firstDerivs, secondDerivs);
      // Evaluate the Module of the strain rate tensor
      sTerm = 0.0;
      for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
        for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
          // Pradtl
          sTerm += (0.5 * (firstDerivs[loopB][loopC] + firstDerivs[loopC][loopB])) * (0.5 * (firstDerivs[loopB][loopC] + firstDerivs[loopC][loopB]));
          // Baldwin and Lomax
          //sTerm += (0.5 * (firstDerivs[loopB][loopC] - firstDerivs[loopC][loopB])) * (0.5 * (firstDerivs[loopB][loopC] - firstDerivs[loopC][loopB]));
        }
      }
      tmp.clear();
      tmp.push_back(density * 0.4187 * currDist * 0.4187 * currDist * sqrt(2.0 * sTerm));
      turbViscosity.push_back(tmp);
      // Increment the counter
      cellCount++;
    }else{
      // Add zero
      tmp.clear();
      tmp.push_back(0.0);
      turbViscosity.push_back(tmp);
    }
  }
}

// ==================================================================
// EXPORT TO POISSON SOLVER ONLY ELEMENTS WITH POSITIVE CONCENTRATION
// ==================================================================
void MRIScan::exportForPoisson(string inputFileName,double density,double viscosity,
                                         MRIThresholdCriteria* threshold,
                                         const MRIDoubleMat& timeDerivs,
                                         bool PPE_IncludeAccelerationTerm,bool PPE_IncludeAdvectionTerm,bool PPE_IncludeDiffusionTerm,bool PPE_IncludeReynoldsTerm,
                                         bool readMuTFromFile, string muTFile, double smagorinskyCoeff){

  printf("\n");
  printf("Acceleration Term included: %s\n",PPE_IncludeAccelerationTerm ? "TRUE":"FALSE");
  printf("Advection Term included:    %s\n",PPE_IncludeAdvectionTerm ? "TRUE":"FALSE");
  printf("Viscous Term included:      %s\n",PPE_IncludeDiffusionTerm ? "TRUE":"FALSE");
  printf("Reynolds Term included:     %s\n",PPE_IncludeReynoldsTerm ? "TRUE":"FALSE");

  // Declare
  FILE* outFile;
  outFile = fopen(inputFileName.c_str(),"w");
  int totAuxNodes = getTotalAuxNodes();
  double qty = 0.0;

  // WRITE THE TOTAL NUMBER OF DOFs FOR THIS PROBLEM
  fprintf(outFile,"NODEDOF %d\n",1);
  fprintf(outFile,"PROBLEM PPE\n");

  // Get the mappings for nodes that need to be used
  MRIIntVec nodeUsageMap;
  MRIIntVec elUsageMap;
  for(int loopA=0;loopA<totAuxNodes;loopA++){
    nodeUsageMap.push_back(-1);
  }

  // Mark Used Nodes
  int currAuxNode = 0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){
      for(int loopB=0;loopB<topology->cellConnections[loopA].size();loopB++){
        currAuxNode = topology->cellConnections[loopA][loopB];
        nodeUsageMap[currAuxNode] = 1;
      }
    }
  }

  // Renumber Nodes in nodeUsageMap
  int usedNodeCount = 0;
  for(int loopA=0;loopA<totAuxNodes;loopA++){
    if(nodeUsageMap[loopA] > 0){
      nodeUsageMap[loopA] = usedNodeCount;
      usedNodeCount++;
    }
  }

  // Build Element Mapping
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    elUsageMap.push_back(-1);
  }
  int elCount = 0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){
      elUsageMap[loopA] = elCount;
      elCount++;
    }
  }

  // ==================
  // SAVE MESH TOPOLOGY
  // ==================
  // SAVE NODE COORDS ONLY FOR NODES NUMBERED IN USEDNODEMAP
  if(topology->totalCells > 0){
    double pos[3];
    for(int loopA=0;loopA<totAuxNodes;loopA++){
      if(nodeUsageMap[loopA] > -1){
        pos[0] = topology->auxNodesCoords[loopA][0];
        pos[1] = topology->auxNodesCoords[loopA][1];
        pos[2] = topology->auxNodesCoords[loopA][2];
        fprintf(outFile,"NODE %d %19.12e %19.12e %19.12e\n",nodeUsageMap[loopA]+1,pos[0],pos[1],pos[2]);
      }
    }

    // ========================
    // SAVE ELEMENT CONNECTIONS
    // ========================
    elCount = 0;
    for(int loopA=0;loopA<topology->totalCells;loopA++){
      qty = cells[loopA].getQuantity(threshold->thresholdQty);
      if(!threshold->MeetsCriteria(qty)){
        fprintf(outFile,"ELEMENT HEXA8 %d 1 ",elCount+1);
        elCount++;
        for(int loopB=0;loopB<topology->cellConnections[loopA].size();loopB++){
          fprintf(outFile,"%d ",nodeUsageMap[topology->cellConnections[loopA][loopB]] + 1);
        }
        fprintf(outFile,"\n");
      }
    }
  }

  // ========================
  // SAVE ELEMENT DIFFUSIVITY
  // ========================
  elCount = 0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){
      fprintf(outFile,"ELDIFF %d ",elCount+1);
      fprintf(outFile,"%e ",1.0);
      fprintf(outFile,"%e ",1.0);
      fprintf(outFile,"%e ",1.0);
      fprintf(outFile,"\n");
      elCount++;
    }
  }  

  // ========================================
  // TEMP: READ TURBULENT VISCOSITY FROM FILE
  // ========================================
  MRIDoubleMat turbViscosity;
  MRIDoubleVec tmp;

  if(PPE_IncludeReynoldsTerm){
    // If the term needs to read by a file
    if(readMuTFromFile){
      printf("Reading Turbulent Viscosity From File...\n");
      MRIUtils::readTableFromFile(muTFile,turbViscosity,false);
      printf("Turbulent Viscosity Points: %d\n",(int)turbViscosity.size());
    }else{
      printf("Using Smagorinsky Lilly Subgrid scale model with coefficient %.4f...\n",smagorinskyCoeff);
      evalSmagorinskyLillyTurbViscosity(density, smagorinskyCoeff, threshold, turbViscosity);
    }
    // Eval the Maximum Module of the turbulent viscosity
    double maxTurbVisc = 0.0;
    for(int loopA=0;loopA<turbViscosity.size();loopA++){
      if(fabs(turbViscosity[loopA][0]) > maxTurbVisc){
        maxTurbVisc = fabs(turbViscosity[loopA][0]);
      }
    }
    printf("Max ABS Turbulent Viscosity: %e\n",maxTurbVisc);
  }else{
     // Zero Turbulent Voscosity
     for(int loopA=0;loopA<topology->totalCells;loopA++){
       tmp.clear();
       tmp.push_back(0.0);
       turbViscosity.push_back(tmp);
     }
  }

  // ==========================================
  // ADD THE TURBULENT VISCOSITY TO THE RESULTS
  // ==========================================
  MRIOutput outMuT("TurbulentViscosity",1);
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Assign To cell
    outMuT.values.push_back(turbViscosity[loopA][0]);
  }
  // Add output quantities
  outputs.push_back(outMuT);

  // ====================================
  // COMPUTE PRESSURE GRADIENT COMPONENTS
  // ====================================
  MRIDoubleMat poissonSourceVec;
  MRIDoubleMat poissonViscousTerm;
  MRIDoubleMat poissonAccelTerm;
  MRIDoubleVec temp;
  MRIDoubleVec tempViscous;
  MRIDoubleVec tempAccel;
  double currValueViscous = 0.0;
  double currValueAdvection = 0.0;
  double currValueAccel = 0.0;

  // ALLOCATE FIRST AND SECOND DERIVATIVES
  MRIDoubleMat firstDerivs;
  MRIDoubleMat secondDerivs;
  firstDerivs.resize(kNumberOfDimensions);
  secondDerivs.resize(kNumberOfDimensions);
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    firstDerivs[loopA].resize(kNumberOfDimensions);
    secondDerivs[loopA].resize(kNumberOfDimensions);
  }

  // ========================================
  // PRINT THE NORM OF THE STRAIN RATE TENSOR
  // ========================================
  MRIOutput outSS("strainRateNorm",1);
  double sTerm = 0.0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Eva Spatial Derivatives
    evalSpaceDerivs(loopA, threshold, firstDerivs, secondDerivs);
    // Evaluate the Module of the velocity Gradient
    sTerm = 0.0;
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
        // Pradtl
        sTerm += (0.5 * (firstDerivs[loopB][loopC] + firstDerivs[loopC][loopB])) * (0.5 * (firstDerivs[loopB][loopC] + firstDerivs[loopC][loopB]));
        // Baldwin and Lomax
        //sTerm += (0.5 * (firstDerivs[loopB][loopC] - firstDerivs[loopC][loopB])) * (0.5 * (firstDerivs[loopB][loopC] - firstDerivs[loopC][loopB]));
      }
    }
    outSS.values.push_back(sqrt(2.0 * sTerm));
  }
  // Add output quantities
  outputs.push_back(outSS);

  // LOOP ON CELLS
  elCount = 0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Only Cells with significant Concentration
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){

      // Eva Spatial Derivatives
      evalSpaceDerivs(loopA, threshold, firstDerivs, secondDerivs);

      // Eval the Convective term for the current Cell
      temp.clear();
      tempViscous.clear();
      tempAccel.clear();
      for(int loopB=0;loopB<kNumberOfDimensions;loopB++){

        // Add Advection term if required
        currValueAdvection = 0.0;
        if(PPE_IncludeAdvectionTerm){
          for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
            currValueAdvection += cells[loopA].velocity[loopC] * firstDerivs[loopC][loopB];
          }
        }

        // Add Viscous term if required
        currValueViscous = 0.0;
        if(PPE_IncludeDiffusionTerm){
          for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
            currValueViscous += secondDerivs[loopC][loopB];
          }
        }

        // Add Acceleration term if required
        if(PPE_IncludeAccelerationTerm){
          currValueAccel = timeDerivs[loopA][loopB];
        }else{
          currValueAccel = 0.0;
        }

        // Add Viscous and Acceleration Components
        temp.push_back(density * currValueAdvection);
        tempViscous.push_back((viscosity + turbViscosity[loopA][0]) * currValueViscous);
        tempAccel.push_back(density * currValueAccel);
      }
      // Add to the global Vector: Only Elements with significant concentration
      poissonSourceVec.push_back(temp);
      poissonViscousTerm.push_back(tempViscous);
      poissonAccelTerm.push_back(tempAccel);
      elCount++;
    }else{
      // Add Zero for the cells with negligible concentration
      temp.clear();
      for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
        temp.push_back(0.0);
      }
      poissonSourceVec.push_back(temp);
      poissonViscousTerm.push_back(temp);
      poissonAccelTerm.push_back(temp);
    }
  }

  // =================================================================================
  // SUM THE CONTRIBUTIONS OF ACCELERATION, ADVECTION, DIFFUSION AND REYNOLDS STRESSES
  // =================================================================================
  MRIDoubleMat termSum;
  double currentValue = 0.0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    temp.clear();    
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      // Set Current Value to Zero
      currentValue = 0.0;
      // All Terms
      if(PPE_IncludeAccelerationTerm){
        currentValue += poissonAccelTerm[loopA][loopB];
      }
      if(PPE_IncludeAdvectionTerm){
        currentValue += poissonSourceVec[loopA][loopB];
      }
      if(PPE_IncludeDiffusionTerm){
        currentValue -= poissonViscousTerm[loopA][loopB];
      }
      // Store Pressure Gradient Components
      temp.push_back(currentValue);
    }
    termSum.push_back(temp);
  }

  // =====================
  // ADD THE RESULT VECTOR
  // ======================
  MRIOutput outF("PressureGradient",3);
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Loop on the dimensions
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      // Assign To cell
      outF.values.push_back(termSum[loopA][loopB]);
    }
  }
  // Add output quantities
  outputs.push_back(outF);

  // ===================
  // FIND FACES ON WALLS
  // ===================
  int* faceCount = new int[topology->faceConnections.size()];
  for(int loopA=0;loopA<topology->faceConnections.size();loopA++){
    faceCount[loopA] = 0;
  }
  bool* isFaceOnWalls = new bool[topology->faceConnections.size()];
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Only Cells with Reasonable Concentration
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){
      for(int loopB=0;loopB<topology->cellFaces[loopA].size();loopB++){
        faceCount[topology->cellFaces[loopA][loopB]]++;
      }
    }
  }
  for(int loopA=0;loopA<topology->faceConnections.size();loopA++){
    if(faceCount[loopA] == 1){
      isFaceOnWalls[loopA] = true;
    }else{
      isFaceOnWalls[loopA] = false;
    }
  }
  delete [] faceCount;

  // Convert Cell Vector to Face Vector
  MRIDoubleVec poissonSourceFaceVec;
  cellToFacePartial(elUsageMap,threshold,termSum,poissonSourceFaceVec);

  // SET TO ZERO THE FACES NOT ON THE BORDER
  //setWallFluxesToZero(isFaceOnWalls,poissonSourceFaceVec);

  // Eval the integral of the divergence over the cell
  MRIDoubleVec cellDivs;
  cellDivs = evalCellDivergences(poissonSourceFaceVec);

  // COMPUTE SOURCE TERMS TO APPLY
  MRIDoubleVec sourcesToApply;
  double currVol = 0.0;
  double SourceSum = 0.0;
  double totalVolume = 0.0;
  double minVolume = std::numeric_limits<double>::max();
  double maxVolume = -std::numeric_limits<double>::max();

  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Evaluate current cell volume
    currVol = evalCellVolume(loopA);    
    // Add Source Term
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){
      SourceSum += cellDivs[loopA];
      totalVolume += currVol;
      if(currVol > maxVolume){
        maxVolume = currVol;
      }
      if(currVol < minVolume){
        minVolume = currVol;
      }
    }
    // Store source term
    sourcesToApply.push_back(cellDivs[loopA]/currVol);
    // TEST
    //sourcesToApply.push_back(cellDivs[loopA]);
  }
  printf("Min Volume: %f, Max Volume: %f\n",minVolume,maxVolume);
  printf("Total Fluid Volume: %f\n",totalVolume);
  printf("Source term summation: %f\n",SourceSum);

  // SAVE ELEMENT SOURCES TO FILE
  elCount = 0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){
      fprintf(outFile,"ELSOURCE %d %19.12e\n",elCount+1,sourcesToApply[loopA]);
      elCount++;
    }
  }

  // ============================
  // CHECK THE DIVERGENCE THEOREM
  // ============================
  int currCell = 0;
  double divSource = 0.0;
  double divNeu = 0.0;
  double sign = 0.0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Evaluate current cell volume
    currVol = evalCellVolume(loopA);
    // Sum Source Contribution
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){
      divSource += cellDivs[loopA];
    }
  }
  double extNormal[3] = {0.0};
  int currFace = 0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){
      for(int loopB=0;loopB<topology->cellFaces[loopA].size();loopB++){
        currFace = topology->cellFaces[loopA][loopB];
        if(isFaceOnWalls[currFace]){
          getExternalFaceNormal(loopA,loopB,extNormal);
          // Get Sign
          sign = 0.0;
          for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
            sign += extNormal[loopB] * topology->faceNormal[currFace][loopB];
          }
          if(sign < 0.0){
            poissonSourceFaceVec[currFace] *= -1.0;
          }
          divNeu += poissonSourceFaceVec[currFace];
        }
      }
    }
  }
  printf("Divergence theorem check: Divergence SUM: %e, Neumann flux SUM: %e\n",divSource,divNeu);

  // =====================
  // SAVE NEUMANN BOUNDARY
  // =====================
  // Loop on the free faces
  for(int loopA=0;loopA<topology->faceCells.size();loopA++){    

    // Check if the face is on the wall
    //if(faceCells[loopA].size() == 1){
    if(isFaceOnWalls[loopA]){

      // Get Current element
      if(topology->faceCells[loopA].size() == 1){
        currCell = topology->faceCells[loopA][0];
      }else{
        qty = cells[topology->faceCells[loopA][0]].getQuantity(threshold->thresholdQty);
        if(!threshold->MeetsCriteria(qty)){
          currCell = topology->faceCells[loopA][0];
        }else{
          currCell = topology->faceCells[loopA][1];
        }
      }

      // Only Cells with significant concentration
      qty = cells[currCell].getQuantity(threshold->thresholdQty);
      if(!threshold->MeetsCriteria(qty)){
        // Print Neumann Condition
        fprintf(outFile,"FACENEUMANN %d ",elUsageMap[currCell] + 1);
        for(int loopB=0;loopB<topology->faceConnections[loopA].size();loopB++){
          fprintf(outFile,"%d ",nodeUsageMap[topology->faceConnections[loopA][loopB]] + 1);
        }
        fprintf(outFile,"%19.12e\n", - poissonSourceFaceVec[loopA]);
      }else{
        printf("PROBLEM!\n");
      }
    }
  }

  // Close Output file
  fclose(outFile);

  // Free Memory
  delete [] isFaceOnWalls;

  printf("\n");
  printf("Poisson Solver File Exported.\n");
}

// ==========================
// CONVERT CELL ARRAY TO FACE
// ==========================
void MRIScan::cellToFace(bool deleteWalls, MRIThresholdCriteria* thresholdCriteria,
                                   MRIDoubleMat cellVec, MRIDoubleVec &faceVec){
  bool   continueToProcess = false;
  double currentValue = 0.0;
  double faceComponent = 0.0;
  int    currentFace = 0;
  double currFaceArea = 0.0;
  bool   checkPassed = false;
  MRIIntVec resID;
  double currVel = 0.0;

  // Get Total Number Of Faces
  int totalFaces = topology->faceConnections.size();

  // Init
  faceVec.resize(totalFaces);
  resID.resize(totalFaces);
  for(int loopA=0;loopA<totalFaces;loopA++){
    faceVec[loopA] = 0.0;
    resID[loopA] = 0;
  }

  // Loop To Assemble Residual Vector
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Check for BC
    if(deleteWalls){
      currentValue = cells[loopA].getQuantity(thresholdCriteria->thresholdQty);
      continueToProcess = !(thresholdCriteria->MeetsCriteria(currentValue));
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
          currVel = cellVec[loopA][loopC];
          faceComponent += currVel * topology->faceNormal[currentFace][loopC];
        }
        // Assemble
        faceVec[currentFace] = faceVec[currentFace] + currFaceArea * faceComponent;
        resID[currentFace]++;
      }
    }
  }
  // Check Faces
  for(int loopA=0;loopA<totalFaces;loopA++){
    if(deleteWalls) checkPassed = (resID[loopA]>2);
    else checkPassed = (resID[loopA]<1)||(resID[loopA]>2);
    if(checkPassed){
      std::string currentMsgs = "Internal: Wrong Face Connectivity, Face: " + MRIUtils::intToStr(loopA)+ "; Connectivity: " + MRIUtils::intToStr(resID[loopA])+".";
      throw new MRIException(currentMsgs.c_str());
    }
  }

  // Divide By the Number Of Faces
  for(int loopA=0;loopA<totalFaces;loopA++){
    if(resID[loopA]>0) faceVec[loopA] = ((double)faceVec[loopA]/(double)resID[loopA]);
    else faceVec[loopA] = 0.0;
  }
}

// ====================================================================
// CONVERT CELL ARRAY TO FACE - ONLY CELLS WITH A CERTAIN CONCENTRATION
// ====================================================================
void MRIScan::cellToFacePartial(MRIIntVec elUsageMap, MRIThresholdCriteria* threshold,
                                          MRIDoubleMat cellVec, MRIDoubleVec &faceVec){
  double faceComponent = 0.0;
  int    currentFace = 0;
  double currFaceArea = 0.0;
  bool   checkNotPassed = false;
  MRIIntVec resID;
  double currVel = 0.0;

  // Get Total Number Of Faces
  int totalFaces = topology->faceConnections.size();

  // Init
  faceVec.resize(totalFaces);
  resID.resize(totalFaces);
  for(int loopA=0;loopA<totalFaces;loopA++){
    faceVec[loopA] = 0.0;
    resID[loopA] = 0;
  }

  // Loop To Assemble Residual Vector
  double qty = 0.0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Check for BC
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){
      // Loop On Faces
      for(int loopB=0;loopB<k3DNeighbors;loopB++){
        // Get Current Face
        currentFace = topology->cellFaces[loopA][loopB];
        // Get Face Area
        currFaceArea = topology->faceArea[currentFace];
        // Get Normal Velocity
        faceComponent = 0.0;
        for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
          currVel = cellVec[loopA][loopC];
          faceComponent += currVel * topology->faceNormal[currentFace][loopC];
        }
        // Assemble
        faceVec[currentFace] = faceVec[currentFace] + currFaceArea * faceComponent;
        resID[currentFace]++;
      }
    }
  }
  // Check Faces
  for(int loopA=0;loopA<totalFaces;loopA++){
    checkNotPassed = (resID[loopA] > 2);
    if(checkNotPassed){
      std::string currentMsgs = "Internal: Wrong Face Connectivity, Face: " + MRIUtils::intToStr(loopA)+ "; Connectivity: " + MRIUtils::intToStr(resID[loopA])+".";
      throw new MRIException(currentMsgs.c_str());
    }
  }
  // Divide By the Number Of Faces
  for(int loopA=0;loopA<totalFaces;loopA++){
    if(resID[loopA]>0) faceVec[loopA] = ((double)faceVec[loopA]/(double)resID[loopA]);
    else faceVec[loopA] = 0.0;
  }
}

// ====================
// EVALUATE CELL VOLUME
// ====================
double MRIScan::evalCellVolume(int cellNumber){
  // Get Integer Indexes
  int intCoords[3];
  mapIndexToCoords(cellNumber,intCoords);
  return topology->cellLengths[0][intCoords[0]] * topology->cellLengths[1][intCoords[1]] * topology->cellLengths[2][intCoords[2]];
}

// ===========================
// PASS INFO FROM PARENT CLASS
// ===========================
void MRIScan::passScanData(MRICommunicator* comm){
  MRIDoubleVec doubleVec;
  if(comm->currProc == 0){
    // Domain Dimension
    doubleVec.push_back(topology->domainSizeMin[0]);
    doubleVec.push_back(topology->domainSizeMin[1]);
    doubleVec.push_back(topology->domainSizeMin[2]);
    doubleVec.push_back(topology->domainSizeMax[0]);
    doubleVec.push_back(topology->domainSizeMax[1]);
    doubleVec.push_back(topology->domainSizeMax[2]);
    doubleVec.push_back(maxVelModule);
  }
  // Pass Data
  comm->passStdDoubleVector(doubleVec);
  // Copy Scan Data  
  if(comm->currProc != 0){
    topology->domainSizeMin[0] = doubleVec[0];
    topology->domainSizeMin[1] = doubleVec[1];
    topology->domainSizeMin[2] = doubleVec[2];
    topology->domainSizeMax[0] = doubleVec[3];
    topology->domainSizeMax[1] = doubleVec[4];
    topology->domainSizeMax[2] = doubleVec[5];
    maxVelModule = doubleVec[6];
  }

  // Exchange Cell Data
  comm->passCellData(topology->totalCells,cells);
}


// ====================
// DISTRIBUTE SCAN DATA
// ====================
void MRIScan::distributeScanData(MRICommunicator* comm){
  // Pass Scan Data  
  passScanData(comm);
  // Exchange Topology Information
  comm->passStdIntVector(topology->cellTotals);
  comm->passStdDoubleMatrix(topology->cellLengths);
  comm->passStdIntMatrix(topology->cellConnections);
  comm->passStdIntMatrix(topology->cellFaces);
  comm->passStdIntMatrix(topology->faceCells);
  comm->passStdIntMatrix(topology->faceConnections);
  comm->passStdIntMatrix(topology->faceEdges);
  comm->passStdDoubleVector(topology->faceArea);
  comm->passStdDoubleMatrix(topology->faceNormal);
  comm->passStdIntMatrix(topology->edgeConnections);
  comm->passStdIntMatrix(topology->edgeFaces);
}

// =========================================
// CHECK IF THE CELL HAS UNTAGGED NEIGHBOURS
// =========================================
bool MRIScan::hasUntaggedNeighbours(int cell,int* cellTags, bool* isTaggable){
  bool result = false;
  int currentFace = 0;
  int nextCell = 0;
  for(int loopA=0;loopA<topology->cellFaces[cell].size();loopA++){
    currentFace = topology->cellFaces[cell][loopA];
    for(int loopB=0;loopB<topology->faceCells[currentFace].size();loopB++){
      nextCell = topology->faceCells[currentFace][loopB];
      result = result || ((isTaggable[nextCell]) && (cellTags[nextCell] == -1));
    }
  }
  // Return the check
  return result;
}

// ==========================
// TAG CELLS USING NEIGHBOURS
// ==========================
void MRIScan::tagByNeighbour(int tag,int* cellTags, bool* isTaggable,int startingCell){
  // Set Starting Cell
  int currentCell = startingCell;
  cellTags[currentCell] = tag;
  MRIIntVec cellList;
  int currentFace = 0;
  int nextCell = 0;

  bool finished = false;
  while(!finished){
    // Explore Cell by Neighbourhood
    for(int loopA=0;loopA<topology->cellFaces[currentCell].size();loopA++){
      currentFace = topology->cellFaces[currentCell][loopA];
      // Tag Faces
      for(int loopB=0;loopB<topology->faceCells[currentFace].size();loopB++){
        nextCell = topology->faceCells[currentFace][loopB];
        if((isTaggable[nextCell]) && (nextCell != currentCell) && (cellTags[nextCell] == -1) ){
          cellTags[nextCell] = tag;
          if(hasUntaggedNeighbours(nextCell,cellTags,isTaggable)){
            cellList.insert(cellList.begin(),nextCell);
          }
        }
      }
    }
    // Update Current Cell
    if(cellList.size() == 0){
      finished = true;
    }else{
      // Get Cell From the tops
      currentCell = cellList[0];
      // Delete First Element
      cellList.erase (cellList.begin());
    }
  }
}

// ================================
// INTERPOLATE DATA ON THE BOUNDARY
// ================================
void MRIScan::interpolateBoundaryVelocities(MRIThresholdCriteria* threshold){
  throw MRIException("ERROR: Not Implemented.\n");
}

// ==============================================
// PROJECT VELOCITY COMPONENT ALONG NORMAL VECTOR
// ==============================================
void MRIScan::projectCellVelocity(int cell,double* normal){
  double normComponent = 0.0;
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    normComponent += normal[loopA] * cells[cell].velocity[loopA];
  }
  cells[cell].velocity[0] = normal[0] * normComponent;
  cells[cell].velocity[1] = normal[1] * normComponent;
  cells[cell].velocity[2] = normal[2] * normComponent;
}

// ==========================
// FIND CELL ON OPPOSITE SIDE
// ==========================
int MRIScan::getOppositeCell(int cell, double* normal){
  int currFace = 0;
  double normalProd = 0.0;
  for(int loopA=0;loopA<topology->cellFaces[cell].size();loopA++){
    // Get Current Face
    currFace = topology->cellFaces[cell][loopA];
    // Get Product
    normalProd = 0.0;
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      normalProd += normal[loopB] * topology->faceNormal[currFace][loopB];
    }
    // Check the Connectivity of the cells
    if((topology->faceCells[currFace].size() > 1) && fabs(normalProd) > 0.5){
      if(topology->faceCells[currFace][0] == cell){
        return topology->faceCells[currFace][1];
      }else{
        return topology->faceCells[currFace][0];
      }
    }
  }
}

// =====================================
// CLEAN VELOCITY COMPONENTS ON BOUNDARY
// =====================================
void MRIScan::cleanNormalComponentOnBoundary(MRIThresholdCriteria* threshold){
  double currFaceNormal[3] = {0.0};
  int currCell = 0;
  int otherCell = 0;
  double qty = 0.0;
  for(int loopA=0;loopA<topology->faceConnections.size();loopA++){
    if(topology->faceCells[loopA].size() == 1){
      //
      currCell = topology->faceCells[loopA][0];
      qty = cells[currCell].getQuantity(threshold->thresholdQty);
      if(!threshold->MeetsCriteria(qty)){
        // Get Normal
        currFaceNormal[0] = topology->faceNormal[loopA][0];
        currFaceNormal[1] = topology->faceNormal[loopA][1];
        currFaceNormal[2] = topology->faceNormal[loopA][2];
        // Project Velocity of Current Cell
        projectCellVelocity(currCell,currFaceNormal);
        // Get Opposite Cells
        otherCell = getOppositeCell(currCell,currFaceNormal);
        // Project Velocity of Current Cell
        projectCellVelocity(otherCell,currFaceNormal);
      }
    }
  }
}

// ==================================================================
// EXPORT TO POISSON SOLVER ONLY ELEMENTS WITH POSITIVE CONCENTRATION
// ==================================================================
void MRIScan::exportForDistancing(string inputFileName,MRIThresholdCriteria* threshold){

  // Declare
  FILE* outFile;
  outFile = fopen(inputFileName.c_str(),"w");
  int totAuxNodes = getTotalAuxNodes();
  double qty = 0.0;

  // WRITE THE NUMBER OF DOFs FOR THIS PROBLEM
  fprintf(outFile,"NODEDOF %d\n",1);
  fprintf(outFile,"PROBLEM DISTANCE\n");

  // Get the mappings for nodes that need to be used
  MRIIntVec nodeUsageMap;
  MRIIntVec elUsageMap;
  for(int loopA=0;loopA<totAuxNodes;loopA++){
    nodeUsageMap.push_back(-1);
  }

  // Mark Used Nodes
  int currAuxNode = 0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){
      for(int loopB=0;loopB<topology->cellConnections[loopA].size();loopB++){
        currAuxNode = topology->cellConnections[loopA][loopB];
        nodeUsageMap[currAuxNode] = 1;
      }
    }
  }

  // Renumber Nodes in nodeUsageMap
  int usedNodeCount = 0;
  for(int loopA=0;loopA<totAuxNodes;loopA++){
    if(nodeUsageMap[loopA] > 0){
      nodeUsageMap[loopA] = usedNodeCount;
      usedNodeCount++;
    }
  }

  // Build Element Mapping
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    elUsageMap.push_back(-1);
  }
  int elCount = 0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){
      elUsageMap[loopA] = elCount;
      elCount++;
    }
  }

  // ==================
  // SAVE MESH TOPOLOGY
  // ==================
  // SAVE NODE COORDS ONLY FOR NODES NUMBERED IN USEDNODEMAP
  if(topology->totalCells > 0){
    double pos[3];
    for(int loopA=0;loopA<totAuxNodes;loopA++){
      if(nodeUsageMap[loopA] > -1){
        pos[0] = topology->auxNodesCoords[loopA][0];
        pos[1] = topology->auxNodesCoords[loopA][1];
        pos[2] = topology->auxNodesCoords[loopA][2];
        fprintf(outFile,"NODE %d %19.12e %19.12e %19.12e\n",nodeUsageMap[loopA]+1,pos[0],pos[1],pos[2]);
      }
    }

    // ========================
    // SAVE ELEMENT CONNECTIONS
    // ========================
    elCount = 0;
    for(int loopA=0;loopA<topology->totalCells;loopA++){
      qty = cells[loopA].getQuantity(threshold->thresholdQty);
      if(!threshold->MeetsCriteria(qty)){
        fprintf(outFile,"ELEMENT HEXA8 %d 1 ",elCount+1);
        elCount++;
        for(int loopB=0;loopB<topology->cellConnections[loopA].size();loopB++){
          fprintf(outFile,"%d ",nodeUsageMap[topology->cellConnections[loopA][loopB]] + 1);
        }
        fprintf(outFile,"\n");
      }
    }
  }

  // ===================
  // FIND FACES ON WALLS
  // ===================
  int* faceCount = new int[topology->faceConnections.size()];
  for(int loopA=0;loopA<topology->faceConnections.size();loopA++){
    faceCount[loopA] = 0;
  }
  MRIBoolVec isFaceOnWalls(topology->faceConnections.size());
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Only Cells with Reasonable Concentration
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){
      for(int loopB=0;loopB<topology->cellFaces[loopA].size();loopB++){
        faceCount[topology->cellFaces[loopA][loopB]]++;
      }
    }
  }
  for(int loopA=0;loopA<topology->faceConnections.size();loopA++){
    if(faceCount[loopA] == 1){
      isFaceOnWalls[loopA] = true;
    }else{
      isFaceOnWalls[loopA] = false;
    }
  }
  delete [] faceCount;

  // ========================
  // SAVE ELEMENT DIFFUSIVITY
  // ========================
  elCount = 0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){
      fprintf(outFile,"ELDIFF %d ",elCount+1);
      fprintf(outFile,"%e ",1.0);
      fprintf(outFile,"%e ",1.0);
      fprintf(outFile,"%e ",1.0);
      fprintf(outFile,"\n");
      elCount++;
    }
  }

  // ==========================================
  // SAVE DIRICHELET CONDITIONS ON THE BOUNDARY
  // ==========================================
  MRIIntVec diricheletNodes;
  diricheletNodes.resize(topology->totalCells);
  long int currCell = 0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    diricheletNodes[loopA] = 0;
  }
  for(int loopA=0;loopA<topology->faceCells.size();loopA++){

    // Check if the face is on the wall
    if(isFaceOnWalls[loopA]){

      // Get Current element
      if (topology->faceCells[loopA].size() == 1){
        currCell = topology->faceCells[loopA][0];
      }else{
        qty = cells[topology->faceCells[loopA][0]].getQuantity(threshold->thresholdQty);
        if(!threshold->MeetsCriteria(qty)){
          currCell = topology->faceCells[loopA][0];
        }else{
          currCell = topology->faceCells[loopA][1];
        }
      }

      // Only Cells with significant concentration
      qty = cells[currCell].getQuantity(threshold->thresholdQty);
      if(!threshold->MeetsCriteria(qty)){
        // Loop through the Nodes
        for(int loopB=0;loopB<topology->faceConnections[loopA].size();loopB++){
          diricheletNodes[nodeUsageMap[topology->faceConnections[loopA][loopB]]]++;
        }
      }else{
        printf("PROBLEM!\n");
      }
    }
  }

  // WRITE DIRICHELET BOUNDARY CONDITIONS
  elCount = 0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    if(diricheletNodes[loopA] > 0){
      fprintf(outFile,"NODEDIRBC %d %19.12e\n",loopA+1,0.0e0);
    }
  }

  // COMPUTE SOURCE TERMS TO APPLY
  MRIDoubleVec sourcesToApply;
  double currVol = 0.0;

  // APPLY SOURCES
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Evaluate current cell volume
    currVol = evalCellVolume(loopA);
    // Add Source Term
    sourcesToApply.push_back(1.0/currVol);
  }

  // SAVE ELEMENT SOURCES TO FILE
  elCount = 0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    qty = cells[loopA].getQuantity(threshold->thresholdQty);
    if(!threshold->MeetsCriteria(qty)){
      fprintf(outFile,"ELSOURCE %d %19.12e\n",elCount+1,sourcesToApply[loopA]);
      elCount++;
    }
  }
  
  // Close Output file
  fclose(outFile);

  printf("\n");
  printf("Distancing Solver File Exported.\n");
}
