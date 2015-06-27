#include <math.h>
#include <string>
#include <limits>
#include <vector>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>

#include "mriCell.h"
#include "mriTypes.h"
#include "mriStructuredScan.h"
#include "mriThresholdCriteria.h"
#include "mriUtils.h"
#include "mriImagedata.h"
#include "mriConstants.h"
#include "mriVolData.h"
#include "mriException.h"
#include "schMessages.h"

using namespace std;

MRIStructuredScan::MRIStructuredScan(double currentTime):MRIScan(currentTime){
  cellTotals.resize(3);
  cellTotals[0] = 0;
  cellTotals[1] = 0;
  cellTotals[2] = 0;
}

// ==============================================
// READS TOKENIZED STRING AND ASSIGNS PLT OPTIONS
// ==============================================
void assignVTKOptions(std::vector<std::string> tokens, vtkStructuredPointsOptionRecord &vtkOptions){
  for(size_t loopA=0;loopA<tokens.size();loopA++){
    if(tokens[loopA] == "ASCII"){
      vtkOptions.isASCII = true;
      vtkOptions.isDefined[0] = true;
    }else if(tokens[loopA].find("STRUCTURED") != string::npos){
      vtkOptions.isValidDataset = true;
      vtkOptions.isDefined[1] = true;
    }else if(tokens[loopA] == "DIMENSIONS"){
      vtkOptions.dimensions[0] = atoi(tokens[loopA+1].c_str());
      vtkOptions.dimensions[1] = atoi(tokens[loopA+2].c_str());
      vtkOptions.dimensions[2] = atoi(tokens[loopA+3].c_str());
      vtkOptions.isDefined[2] = true;
    }else if(tokens[loopA] == "ORIGIN"){
      vtkOptions.origin[0] = atof(tokens[loopA+1].c_str());
      vtkOptions.origin[1] = atof(tokens[loopA+2].c_str());
      vtkOptions.origin[2] = atof(tokens[loopA+3].c_str());
      vtkOptions.isDefined[3] = true;
    }else if(tokens[loopA] == "SPACING"){
      vtkOptions.spacing[0] = atof(tokens[loopA+1].c_str());
      vtkOptions.spacing[1] = atof(tokens[loopA+2].c_str());
      vtkOptions.spacing[2] = atof(tokens[loopA+3].c_str());
      vtkOptions.isDefined[4] = true;
    }
  }
}

// ================
// COPY CONSTRUCTOR
// ================
MRIStructuredScan::MRIStructuredScan(MRIStructuredScan &copyScan):MRIScan(copyScan){
  // COPY THE CELL TOTALS
  cellLengths.resize(3);
  for(int loopA=0;loopA<3;loopA++){
    cellTotals[loopA] = copyScan.cellTotals[loopA];    
    cellLengths[loopA].resize(copyScan.cellLengths[loopA].size());
    // FILL THE LENGTHS
    for(size_t loopB=0;loopB<copyScan.cellLengths[loopA].size();loopB++){
      cellLengths[loopA][loopB] = copyScan.cellLengths[loopA][loopB];
    }
  }
}

// ===================
// GET TOTAL AUX NODES
// ===================
int MRIStructuredScan::getTotalAuxNodes(){
  return (cellTotals[0] + 1)*(cellTotals[1] + 1)*(cellTotals[2] + 1);
}

// ===============
// GET TOTAL FACES
// ===============
int MRIStructuredScan::GetTotalFaces(){
  return cellTotals[0]*cellTotals[1]*(cellTotals[2] + 1)+
         cellTotals[1]*cellTotals[2]*(cellTotals[0] + 1)+
         cellTotals[2]*cellTotals[0]*(cellTotals[1] + 1);
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
};

void MRIStructuredScan::FillPLTHeader(std::vector<std::string> &pltHeader, bool isFirstFile){
  // Clear Vector
  pltHeader.clear();
  if (isFirstFile){
    pltHeader.push_back("TITLE = ""film_cooling.plt""");
    pltHeader.push_back("VARIABLES = ""X/D""");
    pltHeader.push_back("""Y/D""");
    pltHeader.push_back("""Z/D""");
    pltHeader.push_back("""Conc%""");
    pltHeader.push_back("""Vx/Ubulk""");
    pltHeader.push_back("""Vy/Ubulk""");
    pltHeader.push_back("""Vz/Ubulk""");
    if (hasPressureGradient){
      pltHeader.push_back("""pGradx""");
      pltHeader.push_back("""pGrady""");
      pltHeader.push_back("""pGradz""");   
    }
    if (hasRelativePressure){
      pltHeader.push_back("""relPress""");
    }
    if (hasReynoldsStress){
      pltHeader.push_back("""reyXX""");
      pltHeader.push_back("""reyXY""");
      pltHeader.push_back("""reyXZ""");
      pltHeader.push_back("""reyYY""");
      pltHeader.push_back("""reyYZ""");
      pltHeader.push_back("""reyZZ""");
    }
  }
  pltHeader.push_back("ZONE T=""SubZone""");
  pltHeader.push_back(" STRANDID=0, SOLUTIONTIME="+MRIUtils::FloatToStr(scanTime));
  pltHeader.push_back(" I=35, J=113, K=155, ZONETYPE=Ordered");
  pltHeader.push_back(" DATAPACKING=POINT");
  if ((hasPressureGradient)&&(!hasRelativePressure)&&(!hasReynoldsStress)){
    pltHeader.push_back(" DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )");
  }else if ((hasPressureGradient)&&(hasRelativePressure)&&(!hasReynoldsStress)){
    pltHeader.push_back(" DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )");
  }else if ((hasPressureGradient)&&(hasRelativePressure)&&(hasReynoldsStress)){
    pltHeader.push_back(" DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)");
  }else{
    pltHeader.push_back(" DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )");
  }
};

// returns file size in bytes or -1 if not found.
std::ifstream::pos_type GetFileSize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::in | std::ifstream::binary);
    in.seekg(0, std::ifstream::end);
    return in.tellg(); 
};

// Get Statistic String
std::string MRIStructuredScan::WriteStatistics(){

  std::string myresult = "\n";
  myresult += "--------------------------------\n";
  myresult += "FILE STATISTICS\n";
  myresult += "--------------------------------\n";
  myresult += "Total Number Of Cells Read: "+MRIUtils::IntToStr(totalCellPoints)+"\n";
  myresult += "--------------------------------\n";
  myresult += "Total Number Of Coordinate Cells\n";
  myresult += "X Direction: "+MRIUtils::IntToStr(cellTotals[0])+"\n";
  myresult += "Y Direction: "+MRIUtils::IntToStr(cellTotals[1])+"\n";
  myresult += "Z Direction: "+MRIUtils::IntToStr(cellTotals[2])+"\n";
  myresult += "Cells Lengths\n";
  myresult += "X Direction - MIN: "+MRIUtils::FloatToStr(*min_element(cellLengths[0].begin(),cellLengths[0].end()))+"\n";
  myresult += "X Direction - MAX: "+MRIUtils::FloatToStr(*max_element(cellLengths[0].begin(),cellLengths[0].end()))+"\n";
  myresult += "Y Direction - MIN: "+MRIUtils::FloatToStr(*min_element(cellLengths[1].begin(),cellLengths[1].end()))+"\n";
  myresult += "Y Direction - MAX: "+MRIUtils::FloatToStr(*max_element(cellLengths[1].begin(),cellLengths[1].end()))+"\n";
  myresult += "Z Direction - MIN: "+MRIUtils::FloatToStr(*min_element(cellLengths[2].begin(),cellLengths[2].end()))+"\n";
  myresult += "Z Direction - MAX: "+MRIUtils::FloatToStr(*max_element(cellLengths[2].begin(),cellLengths[2].end()))+"\n";
  myresult += "--------------------------------\n";
  myresult += "Domain Size\n";
  myresult += "Minimum X: "+MRIUtils::FloatToStr(domainSizeMin[0])+"\n";
  myresult += "Maximum X: "+MRIUtils::FloatToStr(domainSizeMax[0])+"\n";
  myresult += "Minimum Y: "+MRIUtils::FloatToStr(domainSizeMin[1])+"\n";
  myresult += "Maximum Y: "+MRIUtils::FloatToStr(domainSizeMax[1])+"\n";
  myresult += "Minimum Z: "+MRIUtils::FloatToStr(domainSizeMin[2])+"\n";
  myresult += "Maximum Z: "+MRIUtils::FloatToStr(domainSizeMax[2])+"\n";
  myresult += "--------------------------------\n";
  myresult += "Maximum Velocity Module: "+MRIUtils::FloatToStr(maxVelModule)+"\n";
  myresult += "--------------------------------\n";
  myresult += "\n";
  // Return String
  return myresult;
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
void MRIStructuredScan::ReorderScan(){
  // Determine The Direct and Inverse Permutations
  WriteSchMessage(std::string("Computing Permutation..."));
  std::vector<int> DirectPerm;
  GetGlobalPermutation(DirectPerm);
  WriteSchMessage(std::string("Done.\n"));
  // Reorder Cells
  WriteSchMessage(std::string("Reordering Cells..."));
  ReorderCells(DirectPerm);
  WriteSchMessage(std::string("Done.\n"));
}

// ====================================
// READS HEADER AND ASSIGNS PLT OPTIONS
// ====================================
void assignPLTOptions(std::vector<std::string> tokens, PLTOptionRecord &pltOptions){
  for(size_t loopA=0;loopA<tokens.size();loopA++){
    if(tokens[loopA] == "I"){
      pltOptions.type = pltUNIFORM;
      pltOptions.i = atoi(tokens[loopA+1].c_str());
    }else if(tokens[loopA] == "J"){
      pltOptions.type = pltUNIFORM;
      pltOptions.j = atoi(tokens[loopA+1].c_str());
    }else if(tokens[loopA] == "K"){
      pltOptions.type = pltUNIFORM;
      pltOptions.k = atoi(tokens[loopA+1].c_str());
    }else if(tokens[loopA] == "N"){
      pltOptions.type = pltSTRUCTURED;
      pltOptions.N = atoi(tokens[loopA+1].c_str());
    }else if(tokens[loopA] == "E"){
      pltOptions.type = pltSTRUCTURED;
      pltOptions.E = atoi(tokens[loopA+1].c_str());
    }else if(tokens[loopA] == "FEBLOCK"){
      pltOptions.type = pltSTRUCTURED;
    }
  }
}

// =======================
// READ SCAN FROM PLT FILE
// =======================
void MRIStructuredScan::ReadPltFile(std::string PltFileName, bool DoReorderCells){
  // Init Line Count
  int lineCount = 0;
  totalCellPoints = 0;

  // Initialize Plt Option Record
  PLTOptionRecord pltOptions;

  // Init Domain Limits
  domainSizeMin[0] =  std::numeric_limits<double>::max();
  domainSizeMin[1] =  std::numeric_limits<double>::max();
  domainSizeMin[2] =  std::numeric_limits<double>::max();
  domainSizeMax[0] = -std::numeric_limits<double>::max();
  domainSizeMax[1] = -std::numeric_limits<double>::max();
  domainSizeMax[2] = -std::numeric_limits<double>::max();

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
  WriteSchMessage(std::string("Computing input file size..."));
  int totalLinesInFile = 0;
  std::string Buffer;
  while (std::getline(PltFile,Buffer)){
    if(!foundheader){
      boost::trim(Buffer);
      boost::split(tokenizedString, Buffer, boost::is_any_of(" ,"), boost::token_compress_on);
      areAllFloats = true;
      assignPLTOptions(tokenizedString, pltOptions);
      //for(size_t loopA=0;loopA<tokenizedString.size();loopA++){
      //  areAllFloats = (areAllFloats && (MRIUtils::isFloat(tokenizedString[loopA])));
      //}
      //foundheader = areAllFloats;
      foundheader = headerCount > 100;
      headerCount++;
    }
    // Increase cell and line number
    totalCellPoints++;
    totalLinesInFile++;
  }
  
  // Done: Computing Input File Size
  WriteSchMessage(std::string("Done.\n"));

  // Reset File
  PltFile.clear();
  PltFile.seekg(0, std::ios::beg);
  MRICell myCellPoint;

  // Skip Comments
  std::string* PltFileHeader = new std::string[headerCount];
  for(int loopA=0;loopA<headerCount;loopA++){
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
  std::vector<double> XCoords;
  std::vector<double> YCoords;
  std::vector<double> ZCoords;
  
  // Read All Lines
  int LocalCount = 0;
  int precentProgress = 0;
  int percentCounted = 0;
  int valueCounter = 0;
  int neededValues = 7;
  double* LocalVal = new double[neededValues];
  double* TempVal = new double[neededValues];
  
  // Reading Input File Message
  WriteSchMessage(std::string("Reading input file...\n"));
  
  while (std::getline(PltFile,Buffer)){
    // Read Line
    lineCount++;
    precentProgress = (int)(((double)lineCount/(double)totalLinesInFile)*100);
    if (((precentProgress % 10) == 0)&&((precentProgress / 10) != percentCounted)){
      percentCounted = (precentProgress / 10);
      WriteSchMessage(std::string("Reading..."+MRIUtils::IntToStr(precentProgress)+"\n"));
    }

    // Tokenize Line
    std::vector<std::string> ResultArray = MRIUtils::ExctractSubStringFromBufferMS(Buffer);
    // Store Local Structure
	try{
      // Set Continue
      Continue = true;
      // Check Ratio between ResultArray.size, valueCounter, neededValues
      if(ResultArray.size()+valueCounter < neededValues){
        // Read the whole Result Array
        for(int loopA=0;loopA<ResultArray.size();loopA++){
          LocalVal[loopA] = atof(ResultArray[loopA].c_str());
        }
        Continue = false;
      }else{
        // Read part of the result array
        for(int loopA=0;loopA<neededValues-valueCounter;loopA++){
          LocalVal[loopA] = atof(ResultArray[loopA].c_str());
        }
        // Put the rest in temporary array
        for(int loopA=0;loopA<ResultArray.size()-(neededValues-valueCounter);loopA++){
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
      std::string outString = "WARNING[*] Error Reading Line: "+MRIUtils::IntToStr(lineCount)+"; Line Skipped.\n";
      printf("%s",outString.c_str());
    }
    if (Continue){
      // Update Limits
      // Min
      if (LocalXCoord<domainSizeMin[0]) domainSizeMin[0] = LocalXCoord;
      if (LocalYCoord<domainSizeMin[1]) domainSizeMin[1] = LocalYCoord;
      if (LocalZCoord<domainSizeMin[2]) domainSizeMin[2] = LocalZCoord;
      // Max
      if (LocalXCoord>domainSizeMax[0]) domainSizeMax[0] = LocalXCoord;
      if (LocalYCoord>domainSizeMax[1]) domainSizeMax[1] = LocalYCoord;
      if (LocalZCoord>domainSizeMax[2]) domainSizeMax[2] = LocalZCoord;

      // Update Max Speeds
      if (CurrentModule>maxVelModule) {
        maxVelModule = CurrentModule;
      }

      // Store Node Coords To Find Grid Size
      MRIUtils::InsertInDoubleList(LocalXCoord,TotalXCoords,XCoords);
      MRIUtils::InsertInDoubleList(LocalYCoord,TotalYCoords,YCoords);
      MRIUtils::InsertInDoubleList(LocalZCoord,TotalZCoords,ZCoords);	  

      // Store Velocity/Concentrations
      LocalCount++;
	  
      // Position
      myCellPoint.position[0] = LocalXCoord;
      myCellPoint.position[1] = LocalYCoord;
      myCellPoint.position[2] = LocalZCoord;
      // Conc
      myCellPoint.concentration = LocalConc;
      // Velocity
      myCellPoint.velocity[0] = LocalXVel;
      myCellPoint.velocity[1] = LocalYVel;
      myCellPoint.velocity[2] = LocalZVel;
	  
      // Add to Vector
      cellPoints.push_back(myCellPoint);

      // Set Continue
      Continue = true;
    }
  }
  delete [] LocalVal;
  delete [] TempVal;

  // Set The Effective Number Of Data Read
  totalCellPoints = LocalCount;

  // Store Total Cells
  cellTotals[0] = TotalXCoords;
  cellTotals[1] = TotalYCoords;
  cellTotals[2] = TotalZCoords;

  // Complete To Full Grid: Set To Zero
  totalCellPoints = TotalXCoords * TotalYCoords * TotalZCoords;

  // Set a Zero mtCellPoint
  myCellPoint.position[0] = 0.0;
  myCellPoint.position[1] = 0.0;
  myCellPoint.position[2] = 0.0;
  myCellPoint.concentration = 0.0;
  myCellPoint.velocity[0] = 0.0;
  myCellPoint.velocity[1] = 0.0;
  myCellPoint.velocity[2] = 0.0;

  // Resize CellPoints
  cellPoints.resize(totalCellPoints,myCellPoint);

  // Resize CellLenghts: UNIFORM CASE
  cellLengths.resize(3);
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    cellLengths[loopA].resize(cellTotals[loopA]);
    if(cellTotals[loopA] == 1){
      cellLengths[loopA][0] = 1.0;
    }else{
      for(int loopB=0;loopB<cellTotals[loopA];loopB++){
        cellLengths[loopA][loopB] = fabs(domainSizeMax[loopA]-domainSizeMin[loopA])/(cellTotals[loopA]-1);
      }
    }
  }

  // Finished Reading File
  WriteSchMessage(std::string("File reading completed.\n"));

  // Close File
  PltFile.close();

  // REORDER CELLS
  if (DoReorderCells){
    ReorderScan();
  }

  // CREATING TOPOLOGY
  WriteSchMessage(std::string("\n"));
  WriteSchMessage(std::string("---------------------\n"));
  WriteSchMessage(std::string(" Creating Topology...\n"));
  WriteSchMessage(std::string("---------------------\n"));
  CreateTopology();

  // WRITE STATISTICS
  std::string CurrentStats = WriteStatistics();
  WriteSchMessage(CurrentStats);

}

// ================
// EXPORT TO LSDYNA
// ================
void MRIStructuredScan::ExportToLSDYNA(std::string LSFileName){
  // Export to LS-DYNA
  printf("Exporting File to LSDyna...");
  FILE* LSFile;
  LSFile = fopen(LSFileName.c_str(),"w");
  // Write Header
  fprintf(LSFile,"*KEYWORD\n");
  fprintf(LSFile,"*NODE\n");
  // Set Scale Factor
  double Scale = (1.0/(maxVelModule))*
                 std::max(fabs(domainSizeMax[0]-domainSizeMin[0]),
                 std::max(fabs(domainSizeMax[1]-domainSizeMin[1]),
                     fabs(domainSizeMax[2]-domainSizeMin[2])))*0.05;
  // Create All Nodes
  double CurrentXCoord,CurrentYCoord,CurrentZCoord;
  int TotalNodes = 0;
	double CurrentModule;
  for(int loopA=0;loopA<totalCellPoints;loopA++)
  {
    CurrentModule = sqrt((cellPoints[loopA].velocity[0]*cellPoints[loopA].velocity[0])+
                         (cellPoints[loopA].velocity[1]*cellPoints[loopA].velocity[1])+
                         (cellPoints[loopA].velocity[2]*cellPoints[loopA].velocity[2]));
    // Export Nodes Based On Original Velocity
    if (CurrentModule>kMathZero)
    {
      // Node 1
      TotalNodes++;
      CurrentXCoord = cellPoints[loopA].position[0]-cellPoints[loopA].velocity[0]*Scale;
      CurrentYCoord = cellPoints[loopA].position[1]-cellPoints[loopA].velocity[1]*Scale;
      CurrentZCoord = cellPoints[loopA].position[2]-cellPoints[loopA].velocity[2]*Scale;
      fprintf(LSFile,"%d,%e,%e,%e,0,0\n",TotalNodes,CurrentXCoord,CurrentYCoord,CurrentZCoord);
      // Node 2
      TotalNodes++;
      CurrentXCoord = cellPoints[loopA].position[0]+cellPoints[loopA].velocity[0]*Scale;
      CurrentYCoord = cellPoints[loopA].position[1]+cellPoints[loopA].velocity[1]*Scale;
      CurrentZCoord = cellPoints[loopA].position[2]+cellPoints[loopA].velocity[2]*Scale;
      fprintf(LSFile,"%d,%e,%e,%e,0,0\n",TotalNodes,CurrentXCoord,CurrentYCoord,CurrentZCoord);
    }
  }
  // WRITE ELEMENTS
  fprintf(LSFile,"*ELEMENT_BEAM\n");
  int Count = 0;
	int Node1,Node2;
  for(int loopA=0;loopA<totalCellPoints;loopA++)
  {
    CurrentModule = sqrt((cellPoints[loopA].velocity[0]*cellPoints[loopA].velocity[0])+
                         (cellPoints[loopA].velocity[1]*cellPoints[loopA].velocity[1])+
                         (cellPoints[loopA].velocity[2]*cellPoints[loopA].velocity[2]));
    if (CurrentModule>kMathZero)
    {
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

// -------------
// Export to CSV
// -------------
void MRIStructuredScan::ExportToCSV(std::string FileName)
{
  printf("Exporting to CSV...");
  std::ofstream OutFile;
  OutFile.open(FileName.c_str());
  // Loop On Cells
  for(int loopA=0;loopA<totalCellPoints;loopA++)
  {
    OutFile << cellPoints[loopA].position[0] << cellPoints[loopA].position[1] << cellPoints[loopA].position[2] <<
               cellPoints[loopA].concentration <<
			         cellPoints[loopA].velocity[0] << cellPoints[loopA].velocity[1] << cellPoints[loopA].velocity[2];
  }
  OutFile.close();
  printf("Done\n");
}

// =================
// Export To TECPLOT
// =================
void MRIStructuredScan::ExportToTECPLOT(std::string FileName, bool isFirstFile)
{
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
  FillPLTHeader(PltFileHeader,isFirstFile);
  // Write Header
	std::string LineString;
	std::string compString = "I";
  for(unsigned int loopA=0;loopA<PltFileHeader.size();loopA++)
  {
    boost::trim(PltFileHeader[loopA]);
    LineString = PltFileHeader[loopA];
    if (LineString.substr(0,1) != compString){
      fprintf(outFile,"%s\n",PltFileHeader[loopA].c_str());
    }else{
      fprintf(outFile," I=%d, J=%d, K=%d, ZONETYPE=Ordered\n",cellTotals[0],cellTotals[1],cellTotals[2]);
    }
  }
  // Loop On Cells
  for(int loopA=0;loopA<totalCellPoints;loopA++)
  {
    if ((hasPressureGradient)&&(!hasRelativePressure)&&(!hasReynoldsStress)){
      fprintf(outFile,"%-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n",
                      cellPoints[loopA].position[0],cellPoints[loopA].position[1],cellPoints[loopA].position[2],
                      cellPoints[loopA].concentration,
                      cellPoints[loopA].velocity[0],cellPoints[loopA].velocity[1],cellPoints[loopA].velocity[2],
                      cellPoints[loopA].pressGrad[0],cellPoints[loopA].pressGrad[1],cellPoints[loopA].pressGrad[2]);
    }else if ((hasPressureGradient)&&(hasRelativePressure)&&(!hasReynoldsStress)){
      fprintf(outFile,"%-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n",
                      cellPoints[loopA].position[0],cellPoints[loopA].position[1],cellPoints[loopA].position[2],
                      cellPoints[loopA].concentration,
                      cellPoints[loopA].velocity[0],cellPoints[loopA].velocity[1],cellPoints[loopA].velocity[2],
                      cellPoints[loopA].pressGrad[0],cellPoints[loopA].pressGrad[1],cellPoints[loopA].pressGrad[2],
                      cellPoints[loopA].relPressure);
    }else if ((hasPressureGradient)&&(hasRelativePressure)&&(hasReynoldsStress)){
      fprintf(outFile,"%-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n",
                      cellPoints[loopA].position[0],cellPoints[loopA].position[1],cellPoints[loopA].position[2],
                      cellPoints[loopA].concentration,
                      cellPoints[loopA].velocity[0],cellPoints[loopA].velocity[1],cellPoints[loopA].velocity[2],
                      cellPoints[loopA].pressGrad[0],cellPoints[loopA].pressGrad[1],cellPoints[loopA].pressGrad[2],
                      cellPoints[loopA].relPressure,
                      cellPoints[loopA].ReStress[0],cellPoints[loopA].ReStress[1],cellPoints[loopA].ReStress[2],
                      cellPoints[loopA].ReStress[3],cellPoints[loopA].ReStress[4],cellPoints[loopA].ReStress[5]);
    }else{
      fprintf(outFile,"%-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n",
                      cellPoints[loopA].position[0],cellPoints[loopA].position[1],cellPoints[loopA].position[2],
                      cellPoints[loopA].concentration,
                      cellPoints[loopA].velocity[0],cellPoints[loopA].velocity[1],cellPoints[loopA].velocity[2]);
    }
  }

  // Close Output file
  fclose(outFile);
  
  // Write Done Message
  WriteSchMessage(std::string("Done\n"));
};

// Read Volume File
int MRIStructuredScan::ReadBinVolFile(std::string FileName, MRIVolData &VolData)
{
  // Open File
  FILE *fp = NULL;
  fp = fopen(FileName.c_str(),"rb");
  if (fp == NULL){
    printf("Error: Failed to read volume  %s\n",FileName.c_str());
    return -1;
  }
  // Read Dimensions
  fread(&VolData.GridX, sizeof(int), 1, fp);
  fread(&VolData.GridY, sizeof(int), 1, fp);
  fread(&VolData.GridZ, sizeof(int), 1, fp);
	
  // Read information on Space, Slice Space and Thickness
  fread(&VolData.SpaceX, sizeof(float), 1, fp);
  fread(&VolData.SpaceY, sizeof(float), 1, fp);
  fread(&VolData.SpaceSlice, sizeof(float), 1, fp);
  fread(&VolData.SpaceThick, sizeof(float), 1, fp);

  // Set the size of the voxel information
  int	size = VolData.GridX*VolData.GridY*VolData.GridZ;

  // Read voxel information
  VolData.Voxels = new short[size];
  int len = fread(VolData.Voxels,sizeof(short),size,fp);
  if (len != size)
  {
    printf("Error: Failed to read enough data !\n");
  }
	
  // Close File
  fclose(fp);

  // Return
  return 0;
};
// ------------------------
// Write Binary Volume File
// ------------------------
void WriteBinVolFile(std::string FileName, MRIVolData VolData)
{
	// Write Binary VOL File
  //printf("Writing VOL Binary File...");
	
	// Assign File	
	FILE *fp;
  fp = fopen(FileName.c_str(), "wb");
	
  // Write Grid Size
	fwrite(&VolData.GridX,1,sizeof(int),fp);
	fwrite(&VolData.GridY,1,sizeof(int),fp);
	fwrite(&VolData.GridZ,1,sizeof(int),fp);

  // Read Header
  fwrite(&VolData.SpaceX,1,sizeof(float),fp);
  fwrite(&VolData.SpaceY,1,sizeof(float),fp);
  fwrite(&VolData.SpaceSlice,1,sizeof(float),fp);
  fwrite(&VolData.SpaceThick,1,sizeof(float),fp);

  // Allocate Quantities
  int DataSize = VolData.GridX * VolData.GridY * VolData.GridZ;
  for(int LoopA=0;LoopA<DataSize;LoopA++) fwrite(&VolData.Voxels[LoopA],sizeof(short),1,fp);

  // Close file
	fclose(fp);
  //printf("Done\n");
};

// ========================
// VALIDATE VOL BINARY DATA
// ========================
bool MRIStructuredScan::ValidateVOLBinData(MRIVolData &VolDataAn, MRIVolData &VolDataX, MRIVolData &VolDataY, MRIVolData &VolDataZ){
  // Init Result
  bool result = true;
  // Check Grid Compatibility
  // GridX
  result = MRIUtils::Compare4Integer4(VolDataAn.GridX,VolDataX.GridX,VolDataY.GridX,VolDataZ.GridX);
  if (!result) return result;
  // GridY
  result = MRIUtils::Compare4Integer4(VolDataAn.GridY,VolDataX.GridY,VolDataY.GridY,VolDataZ.GridY);
  if (!result) return result;
  // GridZ
  result = MRIUtils::Compare4Integer4(VolDataAn.GridZ,VolDataX.GridZ,VolDataY.GridZ,VolDataZ.GridZ);
  if (!result) return result;
  // Space X
  result = MRIUtils::Compare4Single(VolDataAn.SpaceX,VolDataX.SpaceX,VolDataY.SpaceX,VolDataZ.SpaceX);
  if (!result) return result;
  // Space Y
  result = MRIUtils::Compare4Single(VolDataAn.SpaceY,VolDataX.SpaceY,VolDataY.SpaceY,VolDataZ.SpaceY);
  if (!result) return result;
  // Space Slice
  result = MRIUtils::Compare4Single(VolDataAn.SpaceSlice,VolDataX.SpaceSlice,VolDataY.SpaceSlice,VolDataZ.SpaceSlice);
  if (!result) return result;
  // Space Thick
  result = MRIUtils::Compare4Single(VolDataAn.SpaceThick,VolDataX.SpaceThick,VolDataY.SpaceThick,VolDataZ.SpaceThick);
  if (!result) return result;
	// Final Message
  return result;
}

// ----------------------------------------------
// Build The Global Data Structure From VOL Files
// ----------------------------------------------
void MRIStructuredScan::FormGlobadDataFromVOL(MRIVolData &VolDataAn, MRIVolData &VolDataX, MRIVolData &VolDataY, MRIVolData &VolDataZ){

  // Allocate Variables
  int coords[kNumberOfDimensions];
  double pos[kNumberOfDimensions];

  // Cells Totals
  cellTotals[0] = VolDataAn.GridX;
  cellTotals[1] = VolDataAn.GridY;
  cellTotals[2] = VolDataAn.GridZ;

  // Cells Length
  cellLengths.resize(3);
  for(int loopB=0;loopB<cellTotals[0];loopB++){
    cellLengths[0][loopB] = VolDataAn.SpaceX;
  }
  for(int loopB=0;loopB<cellTotals[1];loopB++){
    cellLengths[1][loopB] = VolDataAn.SpaceY;
  }
  for(int loopB=0;loopB<cellTotals[2];loopB++){
    cellLengths[2][loopB] = VolDataAn.SpaceSlice;
  }

  // Velocities And Concentrations for all Measure Points
  totalCellPoints = cellTotals[0] * cellTotals[1] * cellTotals[2];
  maxVelModule = 0.0;
  // Allocate
  //cellPoints.reserve(TotalCellPoints+1);
  double currentModule;
  cellPoints.resize(totalCellPoints);
  for(int LoopA=0;LoopA<totalCellPoints;LoopA++)
  {
    // Concentration
    cellPoints[LoopA].concentration = VolDataAn.Voxels[LoopA];
    // Velocity 
    cellPoints[LoopA].velocity[0] = VolDataX.Voxels[LoopA];
    cellPoints[LoopA].velocity[1] = VolDataY.Voxels[LoopA];
    cellPoints[LoopA].velocity[2] = VolDataZ.Voxels[LoopA];
    // Check Max Module
    currentModule = sqrt((cellPoints[LoopA].velocity[0]*cellPoints[LoopA].velocity[0])+
                         (cellPoints[LoopA].velocity[1]*cellPoints[LoopA].velocity[1])+
                         (cellPoints[LoopA].velocity[2]*cellPoints[LoopA].velocity[2]));
    // Get Max Module
    if (currentModule>maxVelModule){
      maxVelModule = currentModule;
    }
  }

  // Init Domain Limits
  domainSizeMin[0] = std::numeric_limits<double>::max();
  domainSizeMin[1] = std::numeric_limits<double>::max();
  domainSizeMin[2] = std::numeric_limits<double>::max();
  domainSizeMax[0] = -std::numeric_limits<double>::max();
  domainSizeMax[1] = -std::numeric_limits<double>::max();
  domainSizeMax[2] = -std::numeric_limits<double>::max();

  // Get The Position From The Index
  for(int LoopA=0;LoopA<totalCellPoints;LoopA++){

    // Map Index To Coords
    MapIndexToCoords(LoopA,coords);

    // Map Integer Coords to Double Coords
    MapCoordsToPosition(coords,false,pos);

    // Store Position
    cellPoints[LoopA].position[0] = pos[0];
    cellPoints[LoopA].position[1] = pos[1];
    cellPoints[LoopA].position[2] = pos[2];

    // Min
    if (pos[0]<domainSizeMin[0]) domainSizeMin[0] = pos[0];
    if (pos[1]<domainSizeMin[1]) domainSizeMin[1] = pos[1];
    if (pos[2]<domainSizeMin[2]) domainSizeMin[2] = pos[2];

    // Max
    if (pos[0]>domainSizeMax[0]) domainSizeMax[0] = pos[0];
    if (pos[1]>domainSizeMax[1]) domainSizeMax[1] = pos[1];
    if (pos[2]>domainSizeMax[2]) domainSizeMax[2] = pos[2];
  }

  // Write Statistics
  std::string Stats = WriteStatistics();
}

// ======================
// CREATE VOL DATA RECORD
// ======================
void MRIStructuredScan::CreateVolDataRecord(int volDataType, MRIVolData &VolData){
  // Cells Totals
  VolData.GridX = cellTotals[0];
  VolData.GridY = cellTotals[1];
  VolData.GridZ = cellTotals[2];

  // Only Uniform Case
  if(!isUniform()){
    throw MRIException("Error. Mesh is not uniform.\n");
  }

  // Cells Length: ASSUME
  VolData.SpaceX = cellLengths[0][0];
  VolData.SpaceY = cellLengths[1][0];
  VolData.SpaceSlice = cellLengths[2][0];
  VolData.SpaceThick = VolData.SpaceSlice;

  // Write Values
  for(int LoopA=0;LoopA<totalCellPoints;LoopA++)
  {
    switch (volDataType) 
		{
      case kVolAnatomy:   
			  VolData.Voxels[LoopA] = short(MRIUtils::round(cellPoints[LoopA].concentration));
				break;
      case kVolVelocityX: 
			  VolData.Voxels[LoopA] = short(MRIUtils::round(cellPoints[LoopA].velocity[0]));
				break;
      case kVolVelocityY: 
			  VolData.Voxels[LoopA] = short(MRIUtils::round(cellPoints[LoopA].velocity[1]));
				break;
      case kVolVelocityZ: 
			  VolData.Voxels[LoopA] = short(MRIUtils::round(cellPoints[LoopA].velocity[2]));
				break;
      case kPressGradX: 
			  VolData.Voxels[LoopA] = short(MRIUtils::round(cellPoints[LoopA].pressGrad[0]));
				break;
      case kPressGradY: 
			  VolData.Voxels[LoopA] = short(MRIUtils::round(cellPoints[LoopA].pressGrad[1]));
				break;
      case kPressGradZ: 
			  VolData.Voxels[LoopA] = short(MRIUtils::round(cellPoints[LoopA].pressGrad[2]));
				break;
      case kRelPressure: 
			  VolData.Voxels[LoopA] = short(MRIUtils::round(cellPoints[LoopA].relPressure));
				break;
        
    }
  }
};

// Assign Name to Vol File
std::string AssignVOLFileName(int fileType, std::string fileName){
  // Create myPath
  //boost::filesystem::path myPath(fileName);
  // Extend File Name
  switch(fileType){
    case kVolAnatomy: 
      //return myPath.filename().string()+"_An"+myPath.extension().string();
      return fileName+"_An.vol";
			break;
    case kVolVelocityX: 
			//return myPath.filename().string()+"_xVel"+myPath.extension().string();
      return fileName+"_xVel.vol";
			break;
    case kVolVelocityY: 
			//return myPath.filename().string()+"_yVel"+myPath.extension().string();
      return fileName+"_yVel.vol";
			break;
    case kVolVelocityZ: 
			//return myPath.filename().string()+"_zVel"+myPath.extension().string();
      return fileName+"_zVel.vol";
			break;    
    case kPressGradX: 
			//return myPath.filename().string()+"_zVel"+myPath.extension().string();
      return fileName+"_xPGrad.vol";
			break;    
    case kPressGradY: 
			//return myPath.filename().string()+"_zVel"+myPath.extension().string();
      return fileName+"_yPGrad.vol";
			break;    
    case kPressGradZ: 
			//return myPath.filename().string()+"_zVel"+myPath.extension().string();
      return fileName+"_zPGrad.vol";
			break;    
    case kRelPressure: 
			//return myPath.filename().string()+"_zVel"+myPath.extension().string();
      return fileName+"_RelP.vol";
			break;    
      
  }
  return fileName;
}

// -------------
// Export to VOL
// -------------
void MRIStructuredScan::ExportToVOL(std::string FileName){
  // Write Message
  WriteSchMessage("Exporting to VOL Files...");
  // Create Vol Data Records
  MRIVolData VolData;
  // Allocate
  VolData.Voxels = new short[totalCellPoints];

  // Anatomy
  CreateVolDataRecord(kVolAnatomy,VolData);
  std::string TempFileName = AssignVOLFileName(kVolAnatomy,FileName);
  WriteBinVolFile(TempFileName,VolData);
  // Velocity X
  CreateVolDataRecord(kVolVelocityX,VolData);
  TempFileName = AssignVOLFileName(kVolVelocityX,FileName);
  WriteBinVolFile(TempFileName,VolData);
  // Velocity Y
  CreateVolDataRecord(kVolVelocityY,VolData);
  TempFileName = AssignVOLFileName(kVolVelocityY,FileName);
  WriteBinVolFile(TempFileName,VolData);
  // Velocity Z
  CreateVolDataRecord(kVolVelocityZ,VolData);
  TempFileName = AssignVOLFileName(kVolVelocityZ,FileName);
  WriteBinVolFile(TempFileName,VolData);
  // If Pressure Gradient is available then EXPORT
  if (hasPressureGradient){
    // Pressure Gradient X
    CreateVolDataRecord(kPressGradX,VolData);
    TempFileName = AssignVOLFileName(kPressGradX,FileName);
    WriteBinVolFile(TempFileName,VolData);
    // Pressure Gradient Y
    CreateVolDataRecord(kPressGradY,VolData);
    TempFileName = AssignVOLFileName(kPressGradY,FileName);
    WriteBinVolFile(TempFileName,VolData);    
    // Pressure Gradient Z
    CreateVolDataRecord(kPressGradZ,VolData);
    TempFileName = AssignVOLFileName(kPressGradZ,FileName);
    WriteBinVolFile(TempFileName,VolData);
  }
  // If Relative Pressure is available then EXPORT
  if (hasRelativePressure){
    // Relative Pressure
    CreateVolDataRecord(kRelPressure,VolData);
    TempFileName = AssignVOLFileName(kRelPressure,FileName);
    WriteBinVolFile(TempFileName,VolData);
  }
	// Print when finished
	WriteSchMessage("Done\n");
  // Deallocate
  delete [] VolData.Voxels;
};

// ========================
// Get Local Adjacent Plane
// ========================
void MRIStructuredScan::GetLocalStarFaces(int StarNum, int CellsX, int CellsY, int &BottomFace, int &TopFace, int &LeftFace, int &RightFace)
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
void MRIStructuredScan::FlushToFile(std::string FileName){
  // Open Output File
	FILE* outFile;
	outFile = fopen(FileName.c_str(),"w");
	// Write Header
  fprintf(outFile,"%-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s\n","Cell Number","PosX","PosY","PosZ","Conc","VelX","VelY","VelZ");
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    fprintf(outFile,"%-15d %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n",loopA,
            cellPoints[loopA].position[0],cellPoints[loopA].position[1],cellPoints[loopA].position[2],
            cellPoints[loopA].concentration,
            cellPoints[loopA].velocity[0],cellPoints[loopA].velocity[1],cellPoints[loopA].velocity[2]);
  }
	// Close Output file
	fclose(outFile);
}

// READ SCAN FROM VOL FILE
void MRIStructuredScan::ReadScanFromVOLFiles(std::string fileNameAn, std::string fileNameX, std::string fileNameY, std::string fileNameZ){
  // Init
  MRIVolData volDataAn;
  MRIVolData volDataX;
  MRIVolData volDataY;
  MRIVolData volDataZ;
  bool continueProcess = false;
  // Get Anatomy Data
  WriteSchMessage(std::string("Reading VOL File..."));
  ReadBinVolFile(fileNameAn,volDataAn);
  WriteSchMessage(std::string("An."));
  // Open X Velocity Component  
  ReadBinVolFile(fileNameX,volDataX);
  WriteSchMessage(std::string("VelX."));
  // Open X Velocity Component
  ReadBinVolFile(fileNameY,volDataY);
  WriteSchMessage(std::string("VelY."));
  // Open X Velocity Component
  ReadBinVolFile(fileNameZ,volDataZ);
  WriteSchMessage(std::string("VelZ."));
  // Validate Quantities Read
  continueProcess = ValidateVOLBinData(volDataAn,volDataX,volDataY,volDataZ);
  WriteSchMessage(std::string("Validation."));
  // Form Global Data Structure
  if (continueProcess){
    FormGlobadDataFromVOL(volDataAn,volDataX,volDataY,volDataZ);
  }
  WriteSchMessage(std::string("Done.\n"));
}

// Eval The Central Cell for the Domain
int MRIStructuredScan::EvalCentralCell(){
  return MapCoordsToIndex(cellTotals[0]/2,cellTotals[1]/2,cellTotals[2]/2);
}

// ASSEMBLE ENCODING MATRIX
void MRIStructuredScan::AssembleEncodingMatrix(int &totalRows, int &totalColumns, double** &Mat){
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
  totalColumns = 3*totalCellPoints;

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
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Eval Neighbours
    faceXPlus =  GetAdjacentFace(loopA,kfacePlusX);
    faceXMinus = GetAdjacentFace(loopA,kfaceMinusX);
    faceYPlus =  GetAdjacentFace(loopA,kfacePlusY);
    faceYMinus = GetAdjacentFace(loopA,kfaceMinusY);
    faceZPlus =  GetAdjacentFace(loopA,kfacePlusZ);
    faceZMinus = GetAdjacentFace(loopA,kfaceMinusZ);
    // Increment Counters
    faceConn[faceXPlus]++;
    faceConn[faceXMinus]++;
    faceConn[faceYPlus]++;
    faceConn[faceYMinus]++;
    faceConn[faceZPlus]++;
    faceConn[faceZMinus]++;
  }
  // LOOP THROUGH THE FACES
  for(int loopA=0;loopA<totalCellPoints;loopA++){

    // EVAL NEIGHBOURS
    faceXPlus =  GetAdjacentFace(loopA,kfacePlusX);
    faceXMinus = GetAdjacentFace(loopA,kfaceMinusX);
    faceYPlus =  GetAdjacentFace(loopA,kfacePlusY);
    faceYMinus = GetAdjacentFace(loopA,kfaceMinusY);
    faceZPlus =  GetAdjacentFace(loopA,kfacePlusZ);
    faceZMinus = GetAdjacentFace(loopA,kfaceMinusZ);

    // EVAL AREAS
    evalCellAreas(loopA,Areas);
    faceXArea = Areas[0];
    faceYArea = Areas[1];
    faceZArea = Areas[2];

    // Eval Column Number
    faceXColumn = loopA;
    faceYColumn = totalCellPoints + loopA;
    faceZColumn = 2*totalCellPoints + loopA;
    // Assembling Terms
    // X PLUS
    if(faceConn[faceXPlus] == 1){
      Mat[faceXPlus] [faceXColumn] += edgeFactor * faceXArea;
    }else if(faceConn[faceXPlus] == 2){
      Mat[faceXPlus] [faceXColumn] += intFactor * faceXArea;
    }else{
      throw MRIMeshCompatibilityException("Invalid Face Connectivity");
    }
    // X MINUS
    if(faceConn[faceXMinus] == 1){
      Mat[faceXMinus][faceXColumn] += edgeFactor * faceXArea;
    }else if(faceConn[faceXMinus] == 2){
      Mat[faceXMinus][faceXColumn] += intFactor * faceXArea;
    }else{
      throw MRIMeshCompatibilityException("Invalid Face Connectivity");
    }
    // Y PLUS
    if(faceConn[faceYPlus] == 1){
      Mat[faceYPlus] [faceYColumn] += edgeFactor * faceYArea;
    }else if(faceConn[faceYPlus] == 2){
      Mat[faceYPlus] [faceYColumn] += intFactor * faceYArea;
    }else{
      throw MRIMeshCompatibilityException("Invalid Face Connectivity");
    }
    // Y MINUS
    if(faceConn[faceYMinus] == 1){
      Mat[faceYMinus][faceYColumn] += edgeFactor * faceYArea;
    }else if(faceConn[faceYMinus] == 2){
      Mat[faceYMinus][faceYColumn] += intFactor * faceYArea;
    }else{
      throw MRIMeshCompatibilityException("Invalid Face Connectivity");
    }
    // Z PLUS
    if(faceConn[faceZPlus] == 1){
      Mat[faceZPlus] [faceZColumn] += edgeFactor * faceZArea;
    }else if(faceConn[faceZPlus] == 2){
      Mat[faceZPlus] [faceZColumn] += intFactor * faceZArea;
    }else{
      throw MRIMeshCompatibilityException("Invalid Face Connectivity");
    }
    // Z MINUS
    if(faceConn[faceZMinus] == 1){
      Mat[faceZMinus][faceZColumn] += edgeFactor * faceZArea;
    }else if(faceConn[faceZMinus] == 2){
      Mat[faceZMinus][faceZColumn] += intFactor * faceZArea;
    }else{
      throw MRIMeshCompatibilityException("Invalid Face Connectivity");
    }
  }
}

// Assemble Decoding Matrix
void MRIStructuredScan::AssembleDecodingMatrix(int &totalRows, int &totalColumns, double** &Mat){
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
  totalRows = 3*totalCellPoints;
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
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Get Areas
    evalCellAreas(loopA,Areas);
    faceXArea = Areas[0];
    faceYArea = Areas[1];
    faceZArea = Areas[2];

    // Eval Neighbours
    faceXPlus =  GetAdjacentFace(loopA,kfacePlusX);
    faceXMinus = GetAdjacentFace(loopA,kfaceMinusX);
    faceYPlus =  GetAdjacentFace(loopA,kfacePlusY);
    faceYMinus = GetAdjacentFace(loopA,kfaceMinusY);
    faceZPlus =  GetAdjacentFace(loopA,kfacePlusZ);
    faceZMinus = GetAdjacentFace(loopA,kfaceMinusZ);
    // Eval Column Number
    faceXRow = loopA;
    faceYRow = totalCellPoints + loopA;
    faceZRow = 2*totalCellPoints + loopA;
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
void MRIStructuredScan::AssignRandomComponent(const int direction,stdRndGenerator &generator){
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    switch (direction){
      case kdirX:
        cellPoints[loopA].velocity[0] = generator();
        break;
      case kdirY:
        cellPoints[loopA].velocity[1] = generator();
        break;
      case kdirZ:
        cellPoints[loopA].velocity[2] = generator();
        break;
    }
  }
}

// Export Velocities To File As Row in the order X,Y,Z
void MRIStructuredScan::ExportVelocitiesToFile(std::string fileName, bool append){
  
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
    for(int loopB=0;loopB<totalCellPoints;loopB++){
      fprintf(outFile,"%e ",cellPoints[loopB].velocity[loopA]);  
    }
  }
  fprintf(outFile,"\n");  
  
  // Close Output file
  fclose(outFile);
}

// CROP SCAN
void MRIStructuredScan::Crop(double* limitBox){
  // Count The Number Of Cells Remaining
  int reminingCells = 0;
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    if (MRIUtils::IsPointInsideBox(cellPoints[loopA].position[0],cellPoints[loopA].position[1],cellPoints[loopA].position[2],limitBox)){
      reminingCells++;
    }
  }
  // Allocate New Cellpoints
  std::vector<MRICell> tempCellPoints;
  tempCellPoints.resize(reminingCells);
  // Initialize Limits
  domainSizeMin[0] =  std::numeric_limits<double>::max();
  domainSizeMin[1] =  std::numeric_limits<double>::max();
  domainSizeMin[2] =  std::numeric_limits<double>::max();
  domainSizeMax[0] = -std::numeric_limits<double>::max();
  domainSizeMax[1] = -std::numeric_limits<double>::max();
  domainSizeMax[2] = -std::numeric_limits<double>::max();
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
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    if (MRIUtils::IsPointInsideBox(cellPoints[loopA].position[0],cellPoints[loopA].position[1],cellPoints[loopA].position[2],limitBox)){
      // Position
      tempCellPoints[tempCount].position[0] = cellPoints[loopA].position[0];
      tempCellPoints[tempCount].position[1] = cellPoints[loopA].position[1];
      tempCellPoints[tempCount].position[2] = cellPoints[loopA].position[2];
      // Change the Domain Size
      // Min X
      if (tempCellPoints[tempCount].position[0]<domainSizeMin[0]){
        domainSizeMin[0] = tempCellPoints[tempCount].position[0];
      }
      // Max X
      if (tempCellPoints[tempCount].position[0]>domainSizeMax[0]){
        domainSizeMax[0] = tempCellPoints[tempCount].position[0];
      }
      // Min Y
      if (tempCellPoints[tempCount].position[1]<domainSizeMin[1]){
        domainSizeMin[1] = tempCellPoints[tempCount].position[1];
      }
      // Max Y
      if (tempCellPoints[tempCount].position[1]>domainSizeMax[1]){
        domainSizeMin[1] = tempCellPoints[tempCount].position[1];
      }
      // Min Z
      if (tempCellPoints[tempCount].position[0]<domainSizeMin[2]){
        domainSizeMin[2] = tempCellPoints[tempCount].position[2];
      }
      // Max Z
      if (tempCellPoints[tempCount].position[0]>domainSizeMax[2]){
        domainSizeMax[2] = tempCellPoints[tempCount].position[2];
      }
      // Store Coordinates
      MRIUtils::InsertInDoubleList(tempCellPoints[tempCount].position[0],TotalXCoords,XCoords);
      MRIUtils::InsertInDoubleList(tempCellPoints[tempCount].position[1],TotalYCoords,YCoords);
      MRIUtils::InsertInDoubleList(tempCellPoints[tempCount].position[2],TotalZCoords,ZCoords);	  
      // Concentration
      tempCellPoints[tempCount].concentration = cellPoints[loopA].concentration;
      // Velocity
      tempCellPoints[tempCount].velocity[0] = cellPoints[loopA].velocity[0];
      tempCellPoints[tempCount].velocity[1] = cellPoints[loopA].velocity[1];
      tempCellPoints[tempCount].velocity[2] = cellPoints[loopA].velocity[2];
      // Get Velocity Module
      currVelModule = sqrt(tempCellPoints[tempCount].velocity[0]*tempCellPoints[tempCount].velocity[0]+
                           tempCellPoints[tempCount].velocity[1]*tempCellPoints[tempCount].velocity[1]+
                           tempCellPoints[tempCount].velocity[2]*tempCellPoints[tempCount].velocity[2]);
      // Store Maximum Velocity Module
      if (currVelModule>maxVelModule) maxVelModule = currVelModule;
      // Pressure Gradient
      tempCellPoints[tempCount].pressGrad[0] = cellPoints[loopA].pressGrad[0];
      tempCellPoints[tempCount].pressGrad[1] = cellPoints[loopA].pressGrad[1];
      tempCellPoints[tempCount].pressGrad[2] = cellPoints[loopA].pressGrad[2];
      // Relative Pressures
      tempCellPoints[tempCount].relPressure = cellPoints[loopA].relPressure;
      // Update Counter
      tempCount++;
    }
  }
  // Store Totals in various directions
  cellTotals[0] = TotalXCoords;
  cellTotals[1] = TotalYCoords;
  cellTotals[2] = TotalZCoords;
  // Store total number Of Cell
  totalCellPoints = reminingCells;
  cellPoints.resize(reminingCells);
  // Copy Back
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Position
    cellPoints[loopA].position[0] = tempCellPoints[loopA].position[0];
    cellPoints[loopA].position[1] = tempCellPoints[loopA].position[1];
    cellPoints[loopA].position[2] = tempCellPoints[loopA].position[2];
    // Concentration
    cellPoints[loopA].concentration = tempCellPoints[loopA].concentration;
    // Velocity
    cellPoints[loopA].velocity[0] = tempCellPoints[loopA].velocity[0];
    cellPoints[loopA].velocity[1] = tempCellPoints[loopA].velocity[1];
    cellPoints[loopA].velocity[2] = tempCellPoints[loopA].velocity[2];
    // Filtered Velocities
    cellPoints[loopA].auxVector[0] = 0.0;
    cellPoints[loopA].auxVector[1] = 0.0;
    cellPoints[loopA].auxVector[2] = 0.0;
    // Pressure Gradients
    cellPoints[loopA].pressGrad[0] = tempCellPoints[loopA].pressGrad[0];
    cellPoints[loopA].pressGrad[1] = tempCellPoints[loopA].pressGrad[1];
    cellPoints[loopA].pressGrad[2] = tempCellPoints[loopA].pressGrad[2];
    // Relative Pressure
    cellPoints[loopA].relPressure = tempCellPoints[loopA].relPressure;
  }
}

// SCALE VELOCITIES
void MRIStructuredScan::ScaleVelocities(double factor){
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    cellPoints[loopA].velocity[0] *= factor;
    cellPoints[loopA].velocity[1] *= factor;
    cellPoints[loopA].velocity[2] *= factor;
  }
  maxVelModule *= factor;
}

// SCALE POSITIONS
void MRIStructuredScan::ScalePositions(double factor){
  // SCALE CELL LENGTHS
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    for(size_t loopB=0;loopB<cellLengths[loopA].size();loopB++){
      cellLengths[loopA][loopB] *= factor;
    }
  }
  // SCALE POSITIONS
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    cellPoints[loopA].position[0] = (cellPoints[loopA].position[0] - domainSizeMin[0])*factor;
    cellPoints[loopA].position[1] = (cellPoints[loopA].position[1] - domainSizeMin[1])*factor;
    cellPoints[loopA].position[2] = (cellPoints[loopA].position[2] - domainSizeMin[2])*factor;
  }
  // SCALE DOMAIN DIMENSIONS
  // Max
  domainSizeMax[0] = (domainSizeMax[0] - domainSizeMin[0])*factor;
  domainSizeMax[1] = (domainSizeMax[1] - domainSizeMin[1])*factor;
  domainSizeMax[2] = (domainSizeMax[2] - domainSizeMin[2])*factor;
  // Min
  domainSizeMin[0] = 0.0;
  domainSizeMin[1] = 0.0;
  domainSizeMin[2] = 0.0;  
}

// =================
// Write to VTK File
// =================
void MRIStructuredScan::ExportToVTK(std::string fileName){

  // Declare
  bool printAux = true;
  double currXCoord = 0.0;
  double currYCoord = 0.0;
  double currZCoord = 0.0;

  // Open Output File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  // Write Header
  fprintf(outFile,"# vtk DataFile Version 2.0\n");
  fprintf(outFile,"Grid Point Model\n");
  fprintf(outFile,"ASCII\n");

  // Write Data Set
  fprintf(outFile,"DATASET RECTILINEAR_GRID\n");
  fprintf(outFile,"DIMENSIONS %d %d %d\n",cellTotals[0],cellTotals[1],cellTotals[2]);

  // Export X Coordinates
  fprintf(outFile,"X_COORDINATES %d double\n",(int)cellLengths[0].size());
  currXCoord = domainSizeMin[0];
  for(size_t loopA=1;loopA<cellLengths[0].size();loopA++){
    fprintf(outFile,"%e ",currXCoord);
    currXCoord += 0.5*(cellLengths[0][loopA-1] + cellLengths[0][loopA]);
  }
  fprintf(outFile,"%e\n",currXCoord);

  // Export Y Coordinates
  fprintf(outFile,"Y_COORDINATES %d double\n",(int)cellLengths[1].size());
  currYCoord = domainSizeMin[1];
  for(size_t loopA=1;loopA<cellLengths[1].size();loopA++){
    fprintf(outFile,"%e ",currYCoord);
    currYCoord += 0.5*(cellLengths[1][loopA-1] + cellLengths[1][loopA]);
  }
  fprintf(outFile,"%e\n",currYCoord);

  // Export Z Coordinates
  fprintf(outFile,"Z_COORDINATES %d double\n",(int)cellLengths[2].size());
  currZCoord = domainSizeMin[2];
  for(size_t loopA=1;loopA<cellLengths[2].size();loopA++){
    fprintf(outFile,"%e ",currZCoord);
    currZCoord += 0.5*(cellLengths[2][loopA-1] + cellLengths[2][loopA]);
  }
  fprintf(outFile,"%e\n",currZCoord);

  // Export Point quantities
  fprintf(outFile,"POINT_DATA %d\n",totalCellPoints);

  // Print Scalar Concentration
  fprintf(outFile,"SCALARS Concentration float 1\n");
  fprintf(outFile,"LOOKUP_TABLE default\n");
  // Print Concentrations
  for (int loopA=0;loopA<totalCellPoints;loopA++){
    fprintf(outFile,"%e\n",cellPoints[loopA].concentration);
  }

  //Print velocity
  fprintf(outFile,"VECTORS Velocity float\n");
  // Print velocity components
  for (int loopA=0;loopA<totalCellPoints;loopA++){
    fprintf(outFile,"%e %e %e\n",cellPoints[loopA].velocity[0],cellPoints[loopA].velocity[1],cellPoints[loopA].velocity[2]);
    //fprintf(outFile,"%e %e %e\n",cellPoints[loopA].filteredVel[0],cellPoints[loopA].filteredVel[1],cellPoints[loopA].filteredVel[2]);
  }

  // Print Pressure Gradient
  if (hasPressureGradient){
    fprintf(outFile,"VECTORS PressureGrad float\n");
    // Print pressure Gradient
    for (int loopA=0;loopA<totalCellPoints;loopA++){
      fprintf(outFile,"%e %e %e\n",cellPoints[loopA].pressGrad[0],cellPoints[loopA].pressGrad[1],cellPoints[loopA].pressGrad[2]);
    }
  }

  // Print Relative Pressure
  if (hasRelativePressure){
    fprintf(outFile,"SCALARS RelPressure float 1\n");
    fprintf(outFile,"LOOKUP_TABLE default\n");
    // Print Relative Pressure
    for (int loopA=0;loopA<totalCellPoints;loopA++){
      fprintf(outFile,"%e\n",cellPoints[loopA].relPressure);
    }
  }

  // Print Relative Pressure
  if (hasRelativePressure){
    fprintf(outFile,"SCALARS GradientMonitor float 1\n");
    fprintf(outFile,"LOOKUP_TABLE default\n");
    // Print Relative Pressure
    for (int loopA=0;loopA<totalCellPoints;loopA++){
      fprintf(outFile,"%e\n",cellPoints[loopA].auxVector[0]);
    }
  }

  // Print Reynolds Stresses
  if (hasReynoldsStress){
    fprintf(outFile,"TENSORS ReynoldsStress float\n");
    // Print Reynolds Stress Tensor
    for (int loopA=0;loopA<totalCellPoints;loopA++){
      fprintf(outFile,"%e %e %e\n",cellPoints[loopA].ReStress[0],cellPoints[loopA].ReStress[1],cellPoints[loopA].ReStress[2]);
      fprintf(outFile,"%e %e %e\n",cellPoints[loopA].ReStress[1],cellPoints[loopA].ReStress[3],cellPoints[loopA].ReStress[4]);
      fprintf(outFile,"%e %e %e\n",cellPoints[loopA].ReStress[2],cellPoints[loopA].ReStress[4],cellPoints[loopA].ReStress[5]);
      fprintf(outFile,"\n");
    }
  }

  // Print outputs
  for(int loopA=0;loopA<outputs.size();loopA++){
    // Print Header
    if(outputs[loopA].totComponents == 1){
      fprintf(outFile,string("SCALARS " + outputs[loopA].name + " float 1\n").c_str());
      fprintf(outFile,"LOOKUP_TABLE default\n");
      for (int loopB=0;loopB<totalCellPoints;loopB++){
        fprintf(outFile,"%e\n",outputs[loopA].values[loopB]);
      }
    }else{
      int count = 0;
      fprintf(outFile,string("VECTORS " + outputs[loopA].name + " float\n").c_str());
      for (int loopB=0;loopB<totalCellPoints;loopB++){
        for(int loopC=0;loopC<3;loopC++){
          fprintf(outFile,"%e ",outputs[loopA].values[count]);
          count++;
        }
        fprintf(outFile,"\n");
      }
    }
  }

  // Print Vortex Criteria
  /*if (printAux){
    fprintf(outFile,"VECTORS VortexCrit float\n");
    // Print pressure Gradient
    for (int loopA=0;loopA<totalCellPoints;loopA++){
      fprintf(outFile,"%e %e %e\n",cellPoints[loopA].auxVector[0],cellPoints[loopA].auxVector[1],cellPoints[loopA].auxVector[2]);
    }
  }*/


  // Close File
  fclose(outFile);
}

// Eval total Number of Vortices
int MRIStructuredScan::EvalTotalVortex(){
  // Init
  int total = 0;
  int totalSlices = 0;
  int totalStars = 0;
  // YZ Planes
  totalSlices = cellTotals[0];
  totalStars = (cellTotals[1]+1)*(cellTotals[2]+1);
  total += totalSlices * totalStars;
  // XZ Planes
  totalSlices = cellTotals[1];
  totalStars = (cellTotals[0]+1)*(cellTotals[2]+1);
  total += totalSlices * totalStars;
  // XY Planes
  totalSlices = cellTotals[2];
  totalStars = (cellTotals[0]+1)*(cellTotals[1]+1);
  total += totalSlices * totalStars;
  // Return Value
  return total;
}

// ======================
// Read List of Row Files
// ======================
void MRIStructuredScan::ReadRAWFileSequence(std::string fileListName){

  // Init
  MRIImageData data;
  std::vector<std::string> fileList;

  // Read File List
  MRIUtils::ReadFileList(fileListName,fileList);

  // Read the first File
  std::string currFile = fileList[0];

  // Read
  ReadRawImage(currFile,data);

  // Set Totals
  cellTotals[0] = data.sizeX;
  cellTotals[1] = data.sizeY;
  cellTotals[2] = fileList.size();

  // Set cell Lengths 1.0 Uniform in all directions
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    for(size_t loopB=0;loopB<cellLengths[loopA].size();loopB++){
      cellLengths[loopA][loopB] = 1.0;
    }
  }

  // Set total Points
  totalCellPoints = cellTotals[0] * cellTotals[1] * cellTotals[2];

  // Intialize Cells
  MRICell myCellPoint;
  myCellPoint.position[0] = 0.0;
  myCellPoint.position[1] = 0.0;
  myCellPoint.position[2] = 0.0;
  myCellPoint.concentration = 0.0;
  myCellPoint.velocity[0] = 0.0;
  myCellPoint.velocity[1] = 0.0;
  myCellPoint.velocity[2] = 0.0;
  // Resize: CHECK!!!
  cellPoints.resize(totalCellPoints,myCellPoint);

  // Fill Concentration with First Image Data
  for(int loopA=0;loopA<data.sizeX*data.sizeY;loopA++){
    cellPoints[loopA].concentration = data.rawData[loopA];
  }

  // Intialize how many cells read so far
  int readSoFar = data.sizeX*data.sizeY;
  // Loop through all other images
  for(unsigned int loopA=1;loopA<fileList.size();loopA++){
    currFile = fileList[loopA];
    // Read
    ReadRawImage(currFile,data);
    // Check Consistency
    if((data.sizeX != cellTotals[0])||(data.sizeY != cellTotals[1])){
      throw new MRIImageException("Error: Image Sequence not Consistent!");
    }
    // Fill Concentration with First Image Data
    for(int loopB=readSoFar;loopB<readSoFar + data.sizeX*data.sizeY;loopB++){
      cellPoints[loopB].concentration = data.rawData[loopB];
    }
    // Increment readSoFar
    readSoFar += data.sizeX*data.sizeY;
  }
}

// ====================
// Read Raw File Header
// ====================
void ReadRawFileHeader(int &sizeX,int &sizeY,int &numberOfBytes,FILE *fp){
  // Read Dimensions
  int wordCount = 0;
  int currentByte = 0;
  std::string fileWord = "";
  while(wordCount<4){
    currentByte = 0;
    fileWord = "";
    while((currentByte != 10)&&(currentByte != 32)){
      fread(&currentByte, 1, 1, fp);
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
int MRIStructuredScan::ReadRawImage(std::string FileName, MRIImageData &data){
  // Open File
  FILE *fp = NULL;
  fp = fopen(FileName.c_str(),"rb");
  if (fp == NULL){
    printf("Error: Failed to read volume  %s\n",FileName.c_str());
    return -1;
  }
  // Read Raw File Header
  int numberOfBytes = 0;
  ReadRawFileHeader(data.sizeX,data.sizeY,numberOfBytes,fp);

  // Set the size of the voxel information
  int size = data.sizeX * data.sizeY;

  // Read voxel information
  data.rawData = new short[size];
  int len = fread(data.rawData,1,size,fp);
  if (len != size)
  {
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
void ReadExpansionFile(std::string fileName,int* tot,
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
  ResultArray = MRIUtils::ExctractSubStringFromBufferMS(Buffer);
  tot[0] = atoi(ResultArray[0].c_str());
  tot[1] = atoi(ResultArray[1].c_str());
  tot[2] = atoi(ResultArray[2].c_str());

  printf("TOTALS %d %d %d\n",tot[0],tot[1],tot[2]);

  // GET CELL X LENGTHS
  int lengthCount = 0;
  while(lengthCount<tot[0]){
    std::getline(inFile,Buffer);
    boost::trim(Buffer);
    ResultArray = MRIUtils::ExctractSubStringFromBufferMS(Buffer);
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
    ResultArray = MRIUtils::ExctractSubStringFromBufferMS(Buffer);
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
    ResultArray = MRIUtils::ExctractSubStringFromBufferMS(Buffer);
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
  ResultArray = MRIUtils::ExctractSubStringFromBufferMS(Buffer);
  minlimits[0] = atof(ResultArray[0].c_str());
  minlimits[1] = atof(ResultArray[1].c_str());
  minlimits[2] = atof(ResultArray[2].c_str());

  // GET MAX LIMITS
  lineCount++;
  std::getline(inFile,Buffer);
  ResultArray = MRIUtils::ExctractSubStringFromBufferMS(Buffer);
  maxlimits[0] = atof(ResultArray[0].c_str());
  maxlimits[1] = atof(ResultArray[1].c_str());
  maxlimits[2] = atof(ResultArray[2].c_str());

  // GET EXPANSION COEFFICIENTS
  std::vector<double> tempExpansion;
  while(std::getline(inFile,Buffer)){
    ResultArray = MRIUtils::ExctractSubStringFromBufferMS(Buffer);
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
void MRIStructuredScan::ReadFromExpansionFile(std::string fileName,bool applyThreshold, int thresholdType,double thresholdRatio){

  // Allocate Variables
  int tot[3];
  std::vector<double> lengthX;
  std::vector<double> lengthY;
  std::vector<double> lengthZ;
  double minlimits[3];
  double maxlimits[3];
  MRIExpansion* exp = NULL;

  // Read Quantities From File
  ReadExpansionFile(fileName,tot,lengthX,lengthY,lengthZ,minlimits,maxlimits,exp);

  // SET UP SCAN QUANTITIES
  // CELL TOTALS
  cellTotals[0] = tot[0];
  cellTotals[1] = tot[1];
  cellTotals[2] = tot[2];

  // CELL LENGTHS
  cellLengths.resize(kNumberOfDimensions);
  // X
  for(size_t loopA=0;loopA<lengthX.size();loopA++){
    cellLengths[0].push_back(lengthX[loopA]);
  }
  // Y
  for(size_t loopA=0;loopA<lengthY.size();loopA++){
    cellLengths[1].push_back(lengthY[loopA]);
  }
  // Z
  for(size_t loopA=0;loopA<lengthZ.size();loopA++){
    cellLengths[2].push_back(lengthZ[loopA]);
  }

  // DIMENSIONS
  // MIN
  domainSizeMin[0] = minlimits[0];
  domainSizeMin[1] = minlimits[1];
  domainSizeMin[2] = minlimits[2];
  // MAX
  domainSizeMax[0] = maxlimits[0];
  domainSizeMax[1] = maxlimits[1];
  domainSizeMax[2] = maxlimits[2];
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
  totalCellPoints = cellTotals[0]*cellTotals[1]*cellTotals[2];
  // CREATE NEW CELL
  MRICell newCell;
  // INITIALIZE QUANTITIES TO ZERO
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // ADD IT TO CELL POINTS
    cellPoints.push_back(newCell);
  }

  // INITIALIZE POSITIONS
  int intCoords[3] = {0};
  double Pos[3] = {0.0};
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    MapIndexToCoords(loopA,intCoords);
    MapCoordsToPosition(intCoords,true,Pos);
    cellPoints[loopA].setQuantity(kQtyPositionX,domainSizeMin[0] + Pos[0]);
    cellPoints[loopA].setQuantity(kQtyPositionY,domainSizeMin[1] + Pos[1]);
    cellPoints[loopA].setQuantity(kQtyPositionZ,domainSizeMin[2] + Pos[2]);
  }

  // REORDER MODEL
  ReorderScan();

  // BUILD TOPOLOGY
  CreateTopology();

  // REBUILD SCAN
  //RebuildFromExpansion(expansion,true);
  // No Constant Flux
  RebuildFromExpansion(expansion,false);
}

// ====================
// WRITE EXPANSION FILE
// ====================
void MRIStructuredScan::WriteExpansionFile(std::string fileName){
  // Open Output File
  FILE* fid;
  fid = fopen(fileName.c_str(),"w");

  // WRITE TOTAL CELLS
  fprintf(fid,"%15d %15d %15d\n",cellTotals[0],cellTotals[1],cellTotals[2]);


  // WRITE CELL LENGTHS X
  for(int loopA=0;loopA<cellTotals[0];loopA++){
    fprintf(fid,"%15.6e\n",cellLengths[0][loopA]);
  }

  // WRITE CELL LENGTHS Y
  for(int loopA=0;loopA<cellTotals[1];loopA++){
    fprintf(fid,"%15.6e\n",cellLengths[1][loopA]);
  }

  // WRITE CELL LENGTHS Z
  for(int loopA=0;loopA<cellTotals[2];loopA++){
    fprintf(fid,"%15.6e\n",cellLengths[2][loopA]);
  }

  // MIN DOMAIN SIZE
  fprintf(fid,"%15.6e %15.6e %15.6e\n",domainSizeMin[0],domainSizeMin[1],domainSizeMin[2]);
  // MAX DOMAIN SIZE
  fprintf(fid,"%15.6e %15.6e %15.6e\n",domainSizeMax[0],domainSizeMax[1],domainSizeMax[2]);

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
void MRIStructuredScan::ThresholdQuantity(int qtyID,double threshold){
  // DECLARE
  WriteSchMessage(std::string("Applying threshold...\n"));
  double centerCellValue = 0.0;
  // LOOP ON ALL CELLS
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Get Value in Current Cell
    centerCellValue = cellPoints[loopA].getQuantity(qtyID);
    // Apply Threshold
    if (centerCellValue<threshold){
      cellPoints[loopA].setQuantity(qtyID,0.0);
    }
  }
}

// ====================================================================
// DETERMINE THREE-DIMENSIONAL COMPONENTS OF THE EXPANSION COEFFICIENTS
// ====================================================================
void MRIStructuredScan::EvalSMPVortexCriteria(MRIExpansion* exp){
  // LOOP ON CELLS
  int idx1 = 0;
  int idx2 = 0;
  int idx3 = 0;
  int idx4 = 0;
  MRIOutput out1("SMPVortexCriterion",3);
  double avVortexIndex = 0.0;
  double totalIntensity = 0.0;
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Loop on the dimensions
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      // Determine the idexes for the adjacent vortices
      getNeighborVortexes(loopA,loopB,idx1,idx2,idx3,idx4);
      // Average the values
      avVortexIndex = 0.25*(exp->vortexCoeff[idx1] + exp->vortexCoeff[idx2] + exp->vortexCoeff[idx3] + exp->vortexCoeff[idx4]);
      // Get total Intensity
      totalIntensity += (exp->vortexCoeff[idx1]*exp->vortexCoeff[idx1] +
                         exp->vortexCoeff[idx2]*exp->vortexCoeff[idx2] +
                         exp->vortexCoeff[idx3]*exp->vortexCoeff[idx3] +
                         exp->vortexCoeff[idx4]*exp->vortexCoeff[idx4]);
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
int MRIStructuredScan::getTotalFaces(){
  return cellTotals[0]*cellTotals[1]*(cellTotals[2]+1)+
         cellTotals[1]*cellTotals[2]*(cellTotals[0]+1)+
         cellTotals[2]*cellTotals[0]*(cellTotals[1]+1);
}

// =======================
// GET NEIGHBORS OF A CELL
// =======================
void MRIStructuredScan::GetNeighbourCells(int CurrentCell,std::vector<int> &cellNeighbors){
  int* coords = new int[3];
  cellNeighbors.clear();
  //Get The Coordinates of the Current Cell
  MapIndexToCoords(CurrentCell,coords);
  // Get Neighbor
  // coords[0]
  if ((coords[0]-1)>=0){
    cellNeighbors.push_back(MapCoordsToIndex(coords[0]-1,coords[1],coords[2]));
  }else{
    cellNeighbors.push_back(-1);
  }
  // coords[1]
  if((coords[0]+1)<cellTotals[0]){
    cellNeighbors.push_back(MapCoordsToIndex(coords[0]+1,coords[1],coords[2]));
  }else{
    cellNeighbors.push_back(-1);
  }
  // coords[2]
  if((coords[1]-1)>=0){
    cellNeighbors.push_back(MapCoordsToIndex(coords[0],coords[1]-1,coords[2]));
  }else{
    cellNeighbors.push_back(-1);
  }
  // coords[3]
  if((coords[1]+1)<cellTotals[1]){
    cellNeighbors.push_back(MapCoordsToIndex(coords[0],coords[1]+1,coords[2]));
  }else{
    cellNeighbors.push_back(-1);
  }
  // coords[4]
  if((coords[2]-1)>=0){
    cellNeighbors.push_back(MapCoordsToIndex(coords[0],coords[1],coords[2]-1));
  }else{
    cellNeighbors.push_back(-1);
  }
    // coords[5]
  if((coords[2]+1)<cellTotals[2]){
    cellNeighbors.push_back(MapCoordsToIndex(coords[0],coords[1],coords[2]+1));
  }else{
    cellNeighbors.push_back(-1);
  }
  // SCRAMBLE VECTOR !!!
  //std::random_shuffle(&cellNeighbors[0], &cellNeighbors[5]);
  // DEALLOCATE
  delete [] coords;
}

// Get Face Starting From Cell and Unit Vector
int MRIStructuredScan::GetFacewithCellVector(int CurrentCell, double* UnitVector){
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

// ==================
// GET CELL FACE AREA
// ==================
void MRIStructuredScan::evalCellAreas(int cellNumber,double* Areas){
  int intCoords[3];
  MapIndexToCoords(cellNumber,intCoords);
  // Get the Three Edge Lengths
  double EdgeX = cellLengths[0][intCoords[0]];
  double EdgeY = cellLengths[1][intCoords[1]];
  double EdgeZ = cellLengths[2][intCoords[2]];
  // Write Results
  Areas[0] = EdgeY * EdgeZ;
  Areas[1] = EdgeX * EdgeZ;
  Areas[2] = EdgeX * EdgeY;
}

// =======================
// BUILD GRID CONNECTIVITY
// =======================
void MRIStructuredScan::buildCellConnections(){
  // Allocate connections
  cellConnections.resize(totalCellPoints);
  // Loop through the cells
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Get Coordinates
    int intCoords[3] = {0};
    int totZSilceNodes = 0;
    int zOffset = 0;
    int yOffset = 0;
    int xOffset = 0;
    int node1,node2,node3,node4;
    int node5,node6,node7,node8;
    // Find Integer Coords
    MapIndexToCoords(loopA,intCoords);
    // Get The Nodes
    // Front Nodes in Z
    // Total Nodes in a Z layer
    totZSilceNodes = (cellTotals[0]+1)*(cellTotals[1]+1);
    zOffset = intCoords[2] * totZSilceNodes;
    yOffset = intCoords[1] * (cellTotals[0]+1);
    xOffset = intCoords[0];
    // Add Node 1
    node1 = (zOffset + yOffset + xOffset);
    cellConnections[loopA].push_back(node1);
    // Add Node 2
    node2 = (zOffset + yOffset + xOffset + 1);
    cellConnections[loopA].push_back(node2);
    // Add Node 3
    node3 = (zOffset + yOffset + xOffset + cellTotals[0] + 1);
    cellConnections[loopA].push_back(node3);
    // Add Node 4
    node4 = (zOffset + yOffset + xOffset + cellTotals[0] + 2);
    cellConnections[loopA].push_back(node4);
    // Change zOffset
    zOffset += (cellTotals[0]+1)*(cellTotals[1]+1);
    // Add Node 5
    node5 = (zOffset + yOffset + xOffset);
    cellConnections[loopA].push_back(node5);
    // Add Node 6
    node6 = (zOffset + yOffset + xOffset + 1);
    cellConnections[loopA].push_back(node6);
    // Add Node 7
    node7 = (zOffset + yOffset + xOffset + cellTotals[0] + 1);
    cellConnections[loopA].push_back(node7);
    // Add Node 8
    node8 = (zOffset + yOffset + xOffset + cellTotals[0] + 2);
    cellConnections[loopA].push_back(node8);
  }
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

// ================================
// GET CELL EXTERNAL NORMAL AT FACE
// ================================
void MRIStructuredScan::getExternalFaceNormal(int cellID, int localFaceID, double* extNormal){
  // Get Face Nodes
  std::vector<int> faceIds;
  getFaceConnections(localFaceID,cellConnections[cellID],faceIds);

  int node1Coords[3] = {0};
  int node2Coords[3] = {0};
  int node3Coords[3] = {0};
  // Get the integer coordinates for the first three nodes
  MapIndexToAuxNodeCoords(faceIds[0],node1Coords);
  MapIndexToAuxNodeCoords(faceIds[1],node2Coords);
  MapIndexToAuxNodeCoords(faceIds[2],node3Coords);
  // Get the positions for the first three nodes
  double node1Pos[3] = {0.0};
  double node2Pos[3] = {0.0};
  double node3Pos[3] = {0.0};
  MapAuxCoordsToPosition(node1Coords,node1Pos);
  MapAuxCoordsToPosition(node2Coords,node2Pos);
  MapAuxCoordsToPosition(node3Coords,node3Pos);
  // Get the difference
  double diff1[3] = {0.0};
  double diff2[3] = {0.0};
  for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
    diff1[loopB] = node2Pos[loopB] - node1Pos[loopB];
    diff2[loopB] = node3Pos[loopB] - node2Pos[loopB];
  }
  // Get the normal
  MRIUtils::Do3DExternalProduct(diff1,diff2,extNormal);
  MRIUtils::Normalize3DVector(extNormal);
}

// ====================================
// GET LOCAL EDGE CONNECTIONS FROM FACE
// ====================================
void getEdgeConnections(int EdgeID, std::vector<int> faceConnections, std::vector<int> &edgeIds){
    edgeIds.clear();
    switch(EdgeID){
      case 0:
        edgeIds.push_back(faceConnections[0]);
        edgeIds.push_back(faceConnections[1]);
        break;
      case 1:
        edgeIds.push_back(faceConnections[1]);
        edgeIds.push_back(faceConnections[2]);
        break;
      case 2:
        edgeIds.push_back(faceConnections[2]);
        edgeIds.push_back(faceConnections[3]);
        break;
      case 3:
        edgeIds.push_back(faceConnections[3]);
        edgeIds.push_back(faceConnections[0]);
        break;
    }
}

// =====================
// ADD FACE TO FACE LIST
// =====================
int MRIStructuredScan::addToFaceConnections(std::vector<std::vector<mriFace* >> &AuxFirstNodeFaceList, std::vector<int> faceIds){
  mriFace* newFace;
  // Get first node in connectivity
  int firstConnectivityNode = *min_element(std::begin(faceIds), std::end(faceIds));
  // Try to find with the first node list
  bool found = false;
  size_t count = 0;
  while((!found)&&(count<AuxFirstNodeFaceList[firstConnectivityNode].size())){
    found = MRIUtils::isSameIntVector(faceIds,AuxFirstNodeFaceList[firstConnectivityNode][count]->connections);
    // Update
    if(!found){
      count++;
    }
  }
  if(!found){
    // Add to Face List
    faceConnections.push_back(faceIds);
    // Add to AuxFirstNodeFaceList
    newFace = new mriFace;
    newFace->number = faceConnections.size()-1;
    for(size_t loopA=0;loopA<faceIds.size();loopA++){
      newFace->connections.push_back(faceIds[loopA]);
    }
    AuxFirstNodeFaceList[firstConnectivityNode].push_back(newFace);
    // Return
    return (faceConnections.size()-1);
  }else{
    return AuxFirstNodeFaceList[firstConnectivityNode][count]->number;
  }
}

// =====================
// ADD EDGE TO FACE LIST
// =====================
int MRIStructuredScan::addToEdgeConnections(std::vector<std::vector<mriEdge*>> &AuxFirstNodeEdgeList, std::vector<int> edgeIds){
  mriEdge* newEdge;
  // Get first node in connectivity
  int firstConnectivityNode = *min_element(std::begin(edgeIds), std::end(edgeIds));
  // Find it in First Node List
  bool found = false;
  size_t count = 0;
  while((!found)&&(count<AuxFirstNodeEdgeList[firstConnectivityNode].size())){
    found = MRIUtils::isSameIntVector(edgeIds,AuxFirstNodeEdgeList[firstConnectivityNode][count]->connections);
    // Update
    if(!found){
      count++;
    }
  }
  if(!found){
    // Add to Edge List
    edgeConnections.push_back(edgeIds);
    // Add to AuxFirstNodeEdgeList
    newEdge = new mriEdge;
    newEdge->number = edgeConnections.size()-1;
    for(size_t loopA=0;loopA<edgeIds.size();loopA++){
      newEdge->connections.push_back(edgeIds[loopA]);
    }
    AuxFirstNodeEdgeList[firstConnectivityNode].push_back(newEdge);
    // Return
    return (edgeConnections.size()-1);
  }else{
    return AuxFirstNodeEdgeList[firstConnectivityNode][count]->number;
  }
}

// =======================
// BUILD FACE CONNECTIVITY
// =======================
void MRIStructuredScan::buildFaceConnections(){
  std::vector<int> faceIds;
  std::vector<std::vector<mriFace* >> AuxFirstNodeFaceList;
  int currFace = 0;
  cellFaces.resize(totalCellPoints);
  AuxFirstNodeFaceList.resize(getTotalAuxNodes());
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    for(int loopB=0;loopB<k3DNeighbors;loopB++){
      // Get Face Connections
      getFaceConnections(loopB,cellConnections[loopA],faceIds);
      // Add to Face Connections
      currFace = addToFaceConnections(AuxFirstNodeFaceList,faceIds);
      // Add to Cell Faces
      cellFaces[loopA].push_back(currFace);
    }
  }
}

// =======================
// BUILD EDGE CONNECTIVITY
// =======================
void MRIStructuredScan::buildFaceCells(){
  faceCells.resize(faceConnections.size());
  int currFace = 0;
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    for(int loopB=0;loopB<cellFaces[loopA].size();loopB++){
      currFace = cellFaces[loopA][loopB];
      faceCells[currFace].push_back(loopA);
    }
  }
}

// =======================
// BUILD EDGE CONNECTIVITY
// =======================
void MRIStructuredScan::buildEdgeConnections(){
  std::vector<int> edgeIds;
  std::vector<std::vector<mriEdge*>> AuxFirstNodeEdgeList;
  int currEdge = 0;
  double coeff = 0.0;
  faceEdges.resize(faceConnections.size());
  AuxFirstNodeEdgeList.resize(getTotalAuxNodes());
  for(size_t loopA=0;loopA<faceConnections.size();loopA++){
    for(int loopB=0;loopB<4;loopB++){
      //if((loopA % 500) == 0){
      //  WriteSchMessage(std::string("Count: ") + std::to_string(loopA) + std::string("\n"));
      //}
      // Get Face Connections
      getEdgeConnections(loopB,faceConnections[loopA],edgeIds);
      // Add to Face Connections
      currEdge = addToEdgeConnections(AuxFirstNodeEdgeList,edgeIds);
      // Add to Face Edges
      faceEdges[loopA].push_back(currEdge);
    }
  }
  // Build edgeFaces
  edgeFaces.resize(edgeConnections.size());
  for(size_t loopA=0;loopA<faceConnections.size();loopA++){
    for(size_t loopB=0;loopB<faceEdges[loopA].size();loopB++){
      currEdge = faceEdges[loopA][loopB];
      // Get Coefficient
      coeff = getEdgeFaceVortexCoeff(currEdge,loopA);
      if(coeff>0.0){
        edgeFaces[currEdge].push_back(loopA+1);
      }else{
        edgeFaces[currEdge].push_back(-(loopA+1));
      }
    }
  }
}

// ======================
// GET VORTEX COEFFICIENT
// ======================
double MRIStructuredScan::getEdgeFaceVortexCoeff(int edgeID, int faceID){
  std::vector<double> edgeDirVector(3);
  std::vector<double> edgeFaceVector(3);
  std::vector<double> resVec(3);
  double res = 0.0;
  // Get Vectors
  getEdgeDirection(edgeID,edgeDirVector);
  getEdgeToFaceDirection(edgeID,faceID,edgeFaceVector);
  // Eval Vector Product
  MRIUtils::Do3DExternalProduct(edgeDirVector,edgeFaceVector,resVec);
  MRIUtils::Normalize3DVector(resVec);
  // Get Sign
  res = 0.0;
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    res = res + resVec[loopA] * faceNormal[faceID][loopA];
  }
  return round(res);
}

// ==================
// GET EDGE DIRECTION
// ==================
void MRIStructuredScan::getEdgeDirection(int edgeID, std::vector<double> &edgeDirVector){
  int node1 = 0;
  int node2 = 0;
  double node1Pos[3] = {0.0};
  double node2Pos[3] = {0.0};
  edgeDirVector.resize(3);

  // Get The Two Nodes
  node1 = edgeConnections[edgeID][0];
  node2 = edgeConnections[edgeID][1];

  // Eval Auxiliary Node Coordinates
  getAuxNodeCoordinates(node1,node1Pos);
  getAuxNodeCoordinates(node2,node2Pos);

  // Get the versor
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    edgeDirVector[loopA] = node1Pos[loopA] - node2Pos[loopA];
  }
  MRIUtils::Normalize3DVector(edgeDirVector);
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    edgeDirVector[loopA] = fabs(edgeDirVector[loopA]);
  }

}

// ==========================
// GET EDGE TO FACE DIRECTION
// ==========================
void MRIStructuredScan::getEdgeToFaceDirection(int edgeID, int faceID, std::vector<double> &edgeFaceVector){
  // Declare
  double ec[3] = {0.0};
  double fc[3] = {0.0};

  // Get Edge Center
  getEdgeCenter(edgeID,ec);

  // Get Face Center
  getFaceCenter(faceID,fc);

  // Get the versor
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    edgeFaceVector[loopA] = fc[loopA] - ec[loopA];
  }
  MRIUtils::Normalize3DVector(edgeFaceVector);
}

// ====================================
// GET COORDINATEDS FOR AUXILIARY NODES
// ====================================
void MRIStructuredScan::getAuxNodeCoordinates(int nodeNum, double* pos){
  int intAuxCoords[3] = {0};
  // Map To Integer Coordinates
  MapIndexToAuxNodeCoords(nodeNum,intAuxCoords);
  // Map To Spatial Position
  MapAuxCoordsToPosition(intAuxCoords,pos);
}

// ===============
// GET EDGE CENTER
// ===============
void MRIStructuredScan::getEdgeCenter(int edgeID, double* ec){
  int node1 = 0;
  int node2 = 0;
  double node1Pos[3] = {0.0};
  double node2Pos[3] = {0.0};

  // Get The Two Nodes
  node1 = edgeConnections[edgeID][0];
  node2 = edgeConnections[edgeID][1];

  // Eval Auxiliary Node Coordinates
  getAuxNodeCoordinates(node1,node1Pos);
  getAuxNodeCoordinates(node2,node2Pos);

  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    ec[loopA] = 0.5*(node1Pos[loopA] + node2Pos[loopA]);
  }
}

// ===============
// GET FACE CENTER
// ===============
void MRIStructuredScan::getFaceCenter(int faceID, double* fc){
  int currNode = 0;
  double pos[3] = {0.0};

  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    fc[loopA] = 0.0;
  }

  for(size_t loopA=0;loopA<faceConnections[faceID].size();loopA++){
    currNode = faceConnections[faceID][loopA];
    getAuxNodeCoordinates(currNode,pos);
    for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
      fc[loopA] += pos[loopA];
    }
  }
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    fc[loopA] /= (double)faceConnections[faceID].size();
  }
}

// ========================================================
// EVAL THE AVERGAGE ERROR BETWEEN FILTEREDVEL AND VELOCITY
// ========================================================
void MRIStructuredScan::RecoverGlobalErrorEstimates(double& AvNormError, double& AvAngleError){
  // Init
  AvNormError = 0.0;
  AvAngleError = 0.0;
  double diffVel[3] = {0.0};
  double diffNorm;
  double cosAlpha,currentAlpha;
  double normVel[3] = {0.0};
  double normFilterVel[3] = {0.0};
  // Loop
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Eval Velocity Difference
    diffVel[0] = cellPoints[loopA].velocity[0]-cellPoints[loopA].auxVector[0];
    diffVel[1] = cellPoints[loopA].velocity[1]-cellPoints[loopA].auxVector[1];
    diffVel[2] = cellPoints[loopA].velocity[2]-cellPoints[loopA].auxVector[2];
    diffNorm = MRIUtils::Do3DEucNorm(diffVel);
    AvNormError = AvNormError + diffNorm;
    // Eval Velocity Angle
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      normVel[loopB] = cellPoints[loopA].velocity[loopB];
      normFilterVel[loopB] = cellPoints[loopA].auxVector[loopB];
    }
    // Normalize
    MRIUtils::Normalize3DVector(normVel);
    MRIUtils::Normalize3DVector(normFilterVel);
    cosAlpha = normVel[0]*normFilterVel[0]+
               normVel[1]*normFilterVel[1]+
               normVel[2]*normFilterVel[2];
    if(cosAlpha>1.0) cosAlpha = 1.0;
    if(cosAlpha<-1.0) cosAlpha = -1.0;
    currentAlpha = acos(cosAlpha);
    AvAngleError = AvAngleError + currentAlpha;
  }
  // Eval Average Values
  if (totalCellPoints>0){
    AvNormError = (AvNormError/totalCellPoints);
    AvAngleError = (AvAngleError/totalCellPoints);
  }else{
    AvNormError = -1.0;
    AvAngleError = -1.0;
  }
}

// ======================
// BUILD FACE AREA VECTOR
// ======================
void MRIStructuredScan::buildFaceAreasAndNormals(){
  double prod = 0.0;
  faceArea.resize(faceConnections.size());
  faceNormal.resize(faceConnections.size());
  double currNormal[3] = {0.0};
  for(size_t loopA=0;loopA<faceConnections.size();loopA++){
    int node1Coords[3] = {0};
    int node2Coords[3] = {0};
    int node3Coords[3] = {0};
    // Get the integer coordinates for the first three nodes
    MapIndexToAuxNodeCoords(faceConnections[loopA][0],node1Coords);
    MapIndexToAuxNodeCoords(faceConnections[loopA][1],node2Coords);
    MapIndexToAuxNodeCoords(faceConnections[loopA][2],node3Coords);
    // Get the positions for the first three nodes
    double node1Pos[3] = {0.0};
    double node2Pos[3] = {0.0};
    double node3Pos[3] = {0.0};
    MapAuxCoordsToPosition(node1Coords,node1Pos);
    MapAuxCoordsToPosition(node2Coords,node2Pos);
    MapAuxCoordsToPosition(node3Coords,node3Pos);   
    // Get the difference
    double diff1[3] = {0.0};
    double diff2[3] = {0.0};
    double diff3[3] = {0.0};
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      diff1[loopB] = node2Pos[loopB] - node1Pos[loopB];
      diff2[loopB] = node3Pos[loopB] - node2Pos[loopB];
      diff3[loopB] = node1Pos[loopB] - cellPoints[faceCells[loopA][0]].position[loopB];
    }
    double d1 = MRIUtils::Do3DEucNorm(diff1);
    double d2 = MRIUtils::Do3DEucNorm(diff2);
    // Evaluate Face Area
    faceArea[loopA] = d1 * d2;
    // Get the normal
    MRIUtils::Do3DExternalProduct(diff1,diff2,currNormal);
    MRIUtils::Normalize3DVector(currNormal);
    if(faceCells[loopA].size() == 1){
      prod = 0.0;
      for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
        prod += currNormal[loopB] * diff3[loopB];
      }
      if(prod>0.0){
        printf("FLIPPED!\n");
        currNormal[0] = - currNormal[0];
        currNormal[1] = - currNormal[1];
        currNormal[2] = - currNormal[2];
      }
    }
    faceNormal[loopA].push_back(currNormal[0]);
    faceNormal[loopA].push_back(currNormal[1]);
    faceNormal[loopA].push_back(currNormal[2]);
  }
}

// ===============================
// CREATE STRUCTURED MESH TOPOLOGY
// ===============================
void MRIStructuredScan::CreateTopology(){
  // Build Cell Connections
  WriteSchMessage(std::string("Build Cell Connection...\n"));
  buildCellConnections();
  // Build Face Connections
  WriteSchMessage(std::string("Build Face Connection...\n"));
  buildFaceConnections();
  buildFaceCells();
  // Build Face Area and Face Normal Vector
  WriteSchMessage(std::string("Build Areas and Normals...\n"));
  buildFaceAreasAndNormals();
  // Build Edge Connections
  WriteSchMessage(std::string("Build Edge Connections...\n"));
  buildEdgeConnections();
  WriteSchMessage(std::string("Topology Creation Completed.\n"));
}

// ===================================
// BUILD TOPOLOGY VECTORS FOR PARMETIS
// ===================================
void MRIStructuredScan::buildMetisConnectivities(int *eptr,int *eind){
  int currNode = 0;
  // Allocate the pointer
  eptr = new int(cellConnections.size());

  // Initialize Counter
  int count = 0;
  // Loop through the cells
  for(int loopA=0;loopA<cellConnections.size();loopA++){
    // Assign the pointer
    eptr[loopA] = count;
    for(int loopB=0;loopB<cellConnections[loopA].size();loopB++){
      // Increment Counter
      count++;
    }
  }
  // Initialize the Connectivities
  eind = new int(count);

  // Loop through the cells
  count = 0;
  for(int loopA=0;loopA<cellConnections.size();loopA++){
    for(int loopB=0;loopB<cellConnections[loopA].size();loopB++){
      // Get Current Node Number
      currNode = cellConnections[loopA][loopB];
      // Increment Counter
      count++;
      // Assign Node
      eind[count] = currNode;
    }
  }
}

// ============================================
// Create Grid from VTK Structured Scan Options
// ============================================
void MRIStructuredScan::CreateGridFromVTKStructuredPoints(vtkStructuredPointsOptionRecord opts){
  // Assign cell totals
  cellTotals[0] = opts.dimensions[0];
  cellTotals[1] = opts.dimensions[1];
  cellTotals[2] = opts.dimensions[2];
  // Assign total number of cells
  totalCellPoints = cellTotals[0] * cellTotals[1] * cellTotals[2];
  // Assign Cell spacing
  cellLengths.resize(3);
  cellLengths[0].resize(cellTotals[0]);
  cellLengths[1].resize(cellTotals[1]);
  cellLengths[2].resize(cellTotals[2]);
  for(int loopA=0;loopA<cellTotals[0];loopA++){
    cellLengths[0][loopA] = opts.spacing[0];
  }
  for(int loopA=0;loopA<cellTotals[1];loopA++){
    cellLengths[1][loopA] = opts.spacing[1];
  }
  for(int loopA=0;loopA<cellTotals[2];loopA++){
    cellLengths[2][loopA] = opts.spacing[2];
  }
  // Set domain size
  // Min
  domainSizeMin[0] = opts.origin[0];
  domainSizeMin[1] = opts.origin[1];
  domainSizeMin[2] = opts.origin[2];
  // Max
  domainSizeMax[0] = opts.origin[0] + (opts.dimensions[0]-1) * opts.spacing[0];
  domainSizeMax[1] = opts.origin[1] + (opts.dimensions[1]-1) * opts.spacing[1];
  domainSizeMax[2] = opts.origin[2] + (opts.dimensions[2]-1) * opts.spacing[2];
  // Allocate the cells
  cellPoints.resize(totalCellPoints);
  MRICell currentCell;
  // Fill the position vectors
  double locCoordX = opts.origin[0];
  double locCoordY = opts.origin[1];
  double locCoordZ = opts.origin[2];
  for(int loopA=0;loopA<cellTotals[2];loopA++){
    for(int loopB=0;loopB<cellTotals[1];loopB++){
      for(int loopC=0;loopC<cellTotals[0];loopC++){
        // Set Cell positions
        currentCell.position[0] = locCoordX;
        currentCell.position[1] = locCoordY;
        currentCell.position[2] = locCoordZ;
        cellPoints.push_back(currentCell);
        locCoordX += opts.spacing[0];
      }
      locCoordY += opts.spacing[1];
    }
    locCoordZ += opts.spacing[2];
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
void InitVTKStructuredPointsOptions(vtkStructuredPointsOptionRecord &opts){
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


// ==========================
// READ VTK STRUCTURED POINTS
// ==========================
void MRIStructuredScan::ReadVTKStructuredPoints(std::string vtkFileName, bool DoReorderCells){
  // Init Line Count
  int lineCount = 0;
  totalCellPoints = 0;

  // Assign File
  WriteSchMessage(std::string("\n"));
  WriteSchMessage(std::string("--- READING STRUCTURED POINT FILE\n"));
  WriteSchMessage(std::string("\n"));
  std::ifstream vtkFile;
  WriteSchMessage(std::string("Open File: ") + vtkFileName + std::string("\n"));
  vtkFile.open(vtkFileName.c_str());

  // Create and initialize vtkOption Vector
  vtkStructuredPointsOptionRecord vtkOptions;
  InitVTKStructuredPointsOptions(vtkOptions);

  // Read Through and look for options
  std::vector<std::string> tokenizedString;
  bool foundAllOptions = false;
  int headerCount = 0;
  WriteSchMessage(std::string("Computing input file size..."));
  int totalLinesInFile = 0;
  std::string Buffer;
  while (std::getline(vtkFile,Buffer)){
    if(!foundAllOptions){
      boost::split(tokenizedString, Buffer, boost::is_any_of(" ,"), boost::token_compress_on);
      // Check if you find options
      assignVTKOptions(tokenizedString, vtkOptions);
      // Check if all options were found
      foundAllOptions = true;
      for(size_t loopA=0;loopA<vtkOptions.numDefined;loopA++){
        foundAllOptions = (foundAllOptions && (vtkOptions.isDefined[loopA]));
      }
      headerCount++;
    }
    // Increase cell and line number
    totalCellPoints++;
    totalLinesInFile++;
  }

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
  CreateGridFromVTKStructuredPoints(vtkOptions);
  WriteSchMessage(std::string("Done.\n"));

  // Reset File
  vtkFile.clear();
  vtkFile.seekg(0, std::ios::beg);

  // Skip Comments
  for(int loopA=0;loopA<headerCount;loopA++){
    std::getline(vtkFile,Buffer);
    lineCount++;
  }

  // Read All Lines
  int LocalCount = 0;
  int precentProgress = 0;
  int percentCounted = 0;
  double currModulus = 0.0;
  double vx; double vy; double vz;
  bool keepProcess;

  // Reading Input File Message
  WriteSchMessage(std::string("Reading input file...\n"));

  while (std::getline(vtkFile,Buffer)){

    // Read Line
    lineCount++;

    // Write Progress
    precentProgress = (int)(((double)lineCount/(double)totalLinesInFile)*100);
    if (((precentProgress % 10) == 0)&&((precentProgress / 10) != percentCounted)){
      percentCounted = (precentProgress / 10);
      WriteSchMessage(std::string("Reading..." + to_string(precentProgress) + "\n"));
    }

    // Tokenize Line
    boost::trim(Buffer);
    boost::split(tokenizedString, Buffer, boost::is_any_of(" ,"), boost::token_compress_on);

    // Store Local Structure
    try{
      if(tokenizedString.size()<3){
        throw "error";
      }
      vx = atof(tokenizedString[0].c_str());
      vy = atof(tokenizedString[1].c_str());
      vz = atof(tokenizedString[2].c_str());
      // GO ahead
      keepProcess = true;
    }catch (...){
      std::string outString = "WARNING[*] Error Reading Line: "+to_string(lineCount)+"; Line Skipped.\n";
      printf("%s",outString.c_str());
      keepProcess = false;
    }

    if(keepProcess){
      cellPoints[LocalCount].velocity[0] = vx;
      cellPoints[LocalCount].velocity[1] = vy;
      cellPoints[LocalCount].velocity[2] = vz;
      currModulus = sqrt((cellPoints[LocalCount].velocity[0]*cellPoints[LocalCount].velocity[0]) +
                         (cellPoints[LocalCount].velocity[1]*cellPoints[LocalCount].velocity[1]) +
                         (cellPoints[LocalCount].velocity[2]*cellPoints[LocalCount].velocity[2]));
      if(currModulus > maxVelModule){
        maxVelModule = currModulus;
      }      
      // Update Counter
      LocalCount++;
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

  // CREATING TOPOLOGY
  WriteSchMessage(std::string("\n"));
  WriteSchMessage(std::string("---------------------\n"));
  WriteSchMessage(std::string(" Creating Topology...\n"));
  WriteSchMessage(std::string("---------------------\n"));
  CreateTopology();

  // WRITE STATISTICS
  std::string CurrentStats = WriteStatistics();
  WriteSchMessage(CurrentStats);

}

// ============================
// GET FACE ID FROM CELL NUMBER
// ============================
int MRIStructuredScan::GetCellFaceID(int CellId,int FaceId){
  int count = 0;
  bool found = false;
  while((!found)&&(count<cellFaces[CellId].size())){
    // Check if found
    found = (cellFaces[CellId][count] == FaceId);
    // Update
    if(!found){
      count++;
    }
  }
  if(!found){
    throw MRIException("ERROR: Face not found in MRIStructuredScan::GetCellFaceID.\n");
  }
  // Return value
  return count;
}

// ========================
// EXPORT TO POISSON SOLVER
// ========================
void MRIStructuredScan::ExportForPOISSON(){
  // Declare
  FILE* outFile;

  // Set File Names
  string nodeFileName("poissonNodes.dat");
  string elementFileName("poissonConnections.dat");
  string diffusivityFileName("poissonDiffusivity.dat");
  string sourceFileName("poissonSources.dat");
  string diricheletFileName("poissonDirBC.dat");
  string bcFileName("poissonFluxBC.dat");

  // TOTAL AUX NODES
  int totAuxNodes = (cellTotals[0]+1)*(cellTotals[1]+1)*(cellTotals[2]+1);

  // ==================
  // SAVE MESH TOPOLOGY
  // ==================

  // SAVE NODE COORDS
  if(totalCellPoints>0){
    outFile = fopen(nodeFileName.c_str(),"w");
    double pos[3];
    for(int loopA=0;loopA<totAuxNodes;loopA++){
      getAuxNodeCoordinates(loopA,pos);
      fprintf(outFile,"%d %e %e %e\n",loopA,pos[0],pos[1],pos[2]);
    }
    // Close Output file
    fclose(outFile);

    // SAVE ELEMENT CONNECTIONS
    outFile = fopen(elementFileName.c_str(),"w");
    // Write Header
    for(int loopA=0;loopA<totalCellPoints;loopA++){
      fprintf(outFile,"%d ",loopA);
      for(int loopB=0;loopB<cellConnections[loopA].size();loopB++){
        fprintf(outFile,"%d ",cellConnections[loopA][loopB]);
      }
      fprintf(outFile,"\n");
    }
    // Close Output file
    fclose(outFile);
  }

  // ==============================
  // SAVE ANISOTROPIC CONDUCIBILITY
  // ==============================
  outFile = fopen(diffusivityFileName.c_str(),"w");
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    fprintf(outFile,"%d ",loopA);
    fprintf(outFile,"%e ",cellPoints[loopA].concentration);
    fprintf(outFile,"%e ",cellPoints[loopA].concentration*0.05);
    fprintf(outFile,"%e ",cellPoints[loopA].concentration*0.05);
    fprintf(outFile,"\n");
  }
  fclose(outFile);

  // ================
  // SAVE SOURCE TERM
  // ================
  MRIDoubleMat poissonSourceVec;
  MRIDoubleMat poissonViscousTerm;
  MRIDoubleVec temp;
  MRIDoubleVec tempViscous;
  double currValue = 0.0;
  double currValueViscous = 0.0;
  // First and Second Derivatives
  double** firstDerivs = new double*[kNumberOfDimensions];
  double** secondDerivs = new double*[kNumberOfDimensions];
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    firstDerivs[loopA] = new double[kNumberOfDimensions];
    secondDerivs[loopA] = new double[kNumberOfDimensions];
  }

  // Loop on cells
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Eva Spatial Derivatives
    EvalSpaceDerivs(loopA, firstDerivs, secondDerivs);
    // Eval the Convective term for the current Cell
    temp.clear();
    tempViscous.clear();
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      currValue = 0.0;
      currValueViscous = 0.0;
      for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
        // Same velocity component on columns
        currValue += cellPoints[loopA].velocity[loopC] * firstDerivs[loopC][loopB];
        currValueViscous += secondDerivs[loopC][loopB];

      }
      // Add Component
      temp.push_back(currValue);
      tempViscous.push_back(4.0e-3 * currValueViscous);
    }
    // Add to the global Vector
    poissonSourceVec.push_back(temp);
    poissonViscousTerm.push_back(tempViscous);
  }

  // Convert Cell Vector to Face Vector
  bool deleteWalls = false;
  MRIThresholdCriteria* crit = new MRIThresholdCriteria(kCriterionLessThen,kNoQuantity,0.0);
  MRIDoubleVec poissonSourceFaceVec;
  cellToFace(deleteWalls,crit,poissonSourceVec,poissonSourceFaceVec);


  // Eval the integral of the divergence over the cell
  MRIDoubleVec cellDivs;
  cellDivs = evalCellDivergences(poissonSourceFaceVec);

  // Divide by the volume
  MRIDoubleVec sourcesToApply;
  double currVol = 0.0;
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Evaluate current cell volume
    currVol = evalCellVolume(loopA);
    // Store source term
    sourcesToApply.push_back(-cellDivs[loopA]/currVol);
  }

  // SAVE ELEMENT SOURCES TO FILE
  outFile = fopen(sourceFileName.c_str(),"w");
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    fprintf(outFile,"%d %e\n",loopA,sourcesToApply[loopA]);
  }
  // Close Output file
  fclose(outFile);

  // =====================
  // SAVE NEUMANN BOUNDARY
  // =====================
  int currCell = 0;
  double currVec[3] = {0.0};
  double normComp = 0.0;
  int faceID = 0.0;
  // Open File
  outFile = fopen(bcFileName.c_str(),"w");
  // Loop on the free faces
  for(int loopA=0;loopA<faceCells.size();loopA++){
    // Check if the face is free
    if(faceCells[loopA].size() == 1){
      // Get Current element
      currCell = faceCells[loopA][0];
      // Get Velocity Component
      currVec[0] = poissonSourceVec[currCell][0] + poissonViscousTerm[currCell][0];
      currVec[1] = poissonSourceVec[currCell][1] + poissonViscousTerm[currCell][1];
      currVec[2] = poissonSourceVec[currCell][2] + poissonViscousTerm[currCell][2];
      // Get Normal Component
      normComp = 0.0;
      for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
        normComp += currVec[loopB] * faceNormal[loopA][loopB];
      }
      // Multiply by the face area
      currValue = faceArea[loopA] * normComp;
      // Print Line
      if(fabs(currValue) > 0.0){
        fprintf(outFile,"%d ",currCell);
        for(int loopB=0;loopB<faceConnections[loopA].size();loopB++){
          fprintf(outFile,"%d ",faceConnections[loopA][loopB]);
        }
        fprintf(outFile,"%e\n",currValue);
      }
    }
  }
  // Close File
  fclose(outFile);

  // =================
  // NO DIRICHELET BCs
  // =================
  // SAVE DIRICHELET CONDITIONS
  //outFile = fopen(diricheletFileName.c_str(),"w");
  // Write Header
  //double pos[3];
  //for(int loopA=0;loopA<totAuxNodes;loopA++){
  //  getAuxNodeCoordinates(loopA,pos);
  //  if((fabs(pos[0]-(domainSizeMin[0]-0.5*cellLengths[0][0]))<1.0e-5)||((fabs(pos[0]-(domainSizeMax[0]-0.5*cellLengths[0][0]))<1.0e-5))){
  //    fprintf(outFile,"%d %e\n",loopA,5.0);
  //  }
  //}
  // Close Output file
  //fclose(outFile);

  // Free Memory
  delete crit;
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    delete [] firstDerivs[loopA];
    delete [] secondDerivs[loopA];
  }
  delete [] firstDerivs;
  delete [] secondDerivs;
}

// ==========================
// CONVERT CELL ARRAY TO FACE
// ==========================
void MRIStructuredScan::cellToFace(bool deleteWalls, MRIThresholdCriteria* thresholdCriteria,
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
  int totalFaces = faceConnections.size();

  // Init
  faceVec.resize(totalFaces);
  resID.resize(totalFaces);
  for(int loopA=0;loopA<totalFaces;loopA++){
    faceVec[loopA] = 0.0;
    resID[loopA] = 0;
  }

  // Loop To Assemble Residual Vector
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Check for BC
    if(deleteWalls){
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
          currVel = cellVec[loopA][loopC];
          faceComponent += currVel * faceNormal[currentFace][loopC];
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
      std::string currentMsgs = "Internal: Wrong Face Connectivity, Face: " + MRIUtils::IntToStr(loopA)+ "; Connectivity: " + MRIUtils::IntToStr(resID[loopA])+".";
      throw new MRIMeshCompatibilityException(currentMsgs.c_str());
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
double MRIStructuredScan::evalCellVolume(int cellNumber){
  // Get Integer Indexes
  int intCoords[3];
  MapIndexToCoords(cellNumber,intCoords);
  return cellLengths[0][intCoords[0]] * cellLengths[1][intCoords[1]] * cellLengths[2][intCoords[2]];
}

// ===========================
// PASS INFO FROM PARENT CLASS
// ===========================
void MRIStructuredScan::passScanData(MRICommunicator* comm){
  MRIDoubleVec doubleVec;
  if(comm->currProc == 0){
    // Domain Dimension
    doubleVec.push_back(domainSizeMin[0]);
    doubleVec.push_back(domainSizeMin[1]);
    doubleVec.push_back(domainSizeMin[2]);
    doubleVec.push_back(domainSizeMax[0]);
    doubleVec.push_back(domainSizeMax[1]);
    doubleVec.push_back(domainSizeMax[2]);
    doubleVec.push_back(maxVelModule);
  }
  // Pass Data
  comm->passStdDoubleVector(doubleVec);
  // Copy Scan Data
  if(comm->currProc != 0){
    domainSizeMin[0] = doubleVec[0];
    domainSizeMin[1] = doubleVec[1];
    domainSizeMin[2] = doubleVec[2];
    domainSizeMax[0] = doubleVec[3];
    domainSizeMax[1] = doubleVec[4];
    domainSizeMax[2] = doubleVec[5];
    maxVelModule = doubleVec[6];
  }

  // Exchange Cell Data
  comm->passCellData(totalCellPoints,cellPoints);
}


// ====================
// DISTRIBUTE SCAN DATA
// ====================
void MRIStructuredScan::DistributeScanData(MRICommunicator* comm){
  // Pass Scan Data
  passScanData(comm);
  // Exchange Topology Information
  comm->passStdIntVector(cellTotals);
  printf("CELL TOTALS: %d %d %d, proc %d\n",cellTotals[0],cellTotals[1],cellTotals[2],comm->currProc);
  comm->passStdDoubleMatrix(cellLengths);
  comm->passStdIntMatrix(cellConnections);
  comm->passStdIntMatrix(cellFaces);
  comm->passStdIntMatrix(faceCells);
  comm->passStdIntMatrix(faceConnections);
  comm->passStdIntMatrix(faceEdges);
  comm->passStdDoubleVector(faceArea);
  comm->passStdDoubleMatrix(faceNormal);
  comm->passStdIntMatrix(edgeConnections);
  comm->passStdIntMatrix(edgeFaces);
}


