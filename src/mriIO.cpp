# include "mriIO.h"

using namespace std;

// ====================================
// READS HEADER AND ASSIGNS PLT OPTIONS
// ====================================
void assignPLTOptions(const MRIStringVec& tokens, pltOptionRecord& pltOptions){
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

// ==============================================
// READS TOKENIZED STRING AND ASSIGNS PLT OPTIONS
// ==============================================
void assignVTKOptions(int lineNum, const MRIStringVec& tokens, vtkStructuredPointsOptionRecord &vtkOptions){
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

// ===========================
// INIT VTK STRUCTURED OPTIONS
// ===========================
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

// ===================
// READ EXPANSION FILE
// ===================
void readExpansionFile(string fileName,
                       MRIIntVec& tot,
                       MRIDoubleVec& lengthX,
                       MRIDoubleVec& lengthY,
                       MRIDoubleVec& lengthZ,
                       MRIDoubleVec& minlimits,
                       MRIDoubleVec& maxlimits,
                       MRIExpansion* exp){

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
  exp->fillFromVector(tempExpansion);

  // CLOSE FILE
  inFile.close();
}

