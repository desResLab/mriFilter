# include "mriUtils.h"

using namespace std;

// Handle the possibiliy of writing messages to some Error Log or to Native Client Instances
void writeSchMessage(std::string Msg){
  printf("%s",Msg.c_str());
}

void writeHeader(){
  writeSchMessage(string("\n"));
  writeSchMessage(string("-----------------------------------\n"));
  writeSchMessage(string(" Flow Manipulation Toolkit\n"));
  writeSchMessage(string(" 2016 - Daniele Schiavazzi,Ph.D.\n"));
  writeSchMessage(string(" Release: 0.2\n"));
  writeSchMessage(string("-----------------------------------\n"));
  writeSchMessage(string("\n"));
}

// MRI UTILITIES
namespace MRIUtils{

// ==================
// WRITE PROGRAM HELP
// ==================
void writeProgramHelp(){
   printf("\n");
   printf("usage mpFilterApp [options]\n");
   printf("Options\n");
   printf("-i --input       Set input file name\n");
   printf("-o --output      Set output file name\n");
   printf("-c --command     Set command file name\n");
   printf("--filter         Apply SMP filter\n");
   printf("--bcfilter       Apply SMP BC filter\n");
   printf("--iterationTol   Set Iteration tolerance (default: 1.0e-3)\n");
   printf("--maxIterations  Set maximum number of iteration (default: 2000)\n");
   printf("--thresholdQty   Set threshold quantity\n");
   printf("                 -1 = No Threshold\n");
   printf("                  0 = X Coordinate\n");
   printf("                  1 = Y Coordinate\n");
   printf("                  2 = Z Coordinate\n");
   printf("                  3 = Concentration\n");
   printf("                  4 = X Velocity\n");
   printf("                  5 = Y Velocity\n");
   printf("                  6 = Z Velocity\n");
   printf("--thresholdType   Set threshold type\n");
   printf("                 -1 = Less than\n");
   printf("                  0 = Greater than\n");
   printf("                  1 = ABS Less than\n");
   printf("                  2 = ABS Greater than\n");
   printf("--iformat        Set input file format\n");
   printf("                  0 = TECPLOT File Format\n");
   printf("                  1 = VTK File Format\n");
   printf("--oformat        Set output file format\n");
   printf("                  0 = TECPLOT File Format\n");
   printf("                  1 = VTK File Format\n");
   printf("\n");
}

// =========================
// CONVERT INTEGER TO STRING
// =========================
string intToStr(int number){
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

// ========================
// CONVERT DOUBLE TO STRING
// ========================
string floatToStr(double number){
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

// ===============================================
// EXTRACT STRING FROM BUFFER WITH MULTIPLE SPACES
// ===============================================
MRIStringVec extractSubStringFromBufferMS(std::string Buffer){
  // Tokenize
  char* BufferCH = const_cast<char*> (Buffer.c_str());
  char* tok = strtok (BufferCH," ");
  std::vector<std::string> parts;
  // Parse
  while ( tok != NULL ) 
  {
    parts.push_back( tok );
    tok = strtok (NULL," ");
	//string myString = tok + "\n";
	//printf(myString.c_str());
  }
  return parts;
}

// =======================
// EVAL TWO-NORM OF VECTOR
// =======================
double do3DEucNorm(const MRIDoubleVec& v){
  double norm2 = 0.0;
  for(int i = 0;i<kNumberOfDimensions;i++){
    norm2 += v[i]*v[i];
  }
  return sqrt(norm2);
}


// ===================
// NORMALIZE 3D VECTOR
// ===================
void normalize3DVector(MRIDoubleVec& v){
  double norm = do3DEucNorm(v);
  if(norm>kMathZero){
    for(int LoopA=0;LoopA<kNumberOfDimensions;LoopA++){
      v[LoopA] = v[LoopA]/norm;
    }
  }
}

// =============================
// INSERT IN GENERIC DOUBLE LIST
// =============================
void insertInList(double value, MRIDoubleVec& list){
  double distance;	  
  for(int LoopA=0;LoopA<list.size();LoopA++){
    distance = fabs(list[LoopA] - value);
    if (distance < kMathZero) return;
  }
  // Insert New value
  list.push_back(value);
}

// ==============================
// INSERT IN GENERIC INTEGER LIST
// ==============================
void insertInList(int value, MRIIntVec& list){
  for(int LoopA=0;LoopA<list.size();LoopA++){
    if (list[LoopA] == value) return;
  }
  // Insert New value
  list.push_back(value);
}


// ==========================
// GENERATE STANDARD GAUSSIAN
// ==========================
double generateStandardGaussian(double stDev){
  // Allocate Vector
  boost::random::normal_distribution<> dist(0.0,stDev);
  // Add Random Component
  return dist(realGen);
}

// ============
// ROUND VALUES
// ============
int round(double value){
  return (int)floor(value+0.5);
}

// ==========================
// WRITE VECTOR GRAPH TO FILE
// ==========================
void writeGraphToFile(std::string fileName, int vecSize, std::vector<double> &vecX, std::vector<double> &vecY){
  // Open Output Files
	FILE* outFile;
	outFile = fopen(fileName.c_str(),"w");
	// Write Header
  for(int loopA=0;loopA<vecSize;loopA++){
    fprintf(outFile,"%d %e %e\n",loopA,vecX[loopA],vecY[loopA]);
  }
	// Close Output file
	fclose(outFile);  
}

// ========================
// PERFORM EXTERNAL PRODUCT
// ========================
void do3DExternalProduct(const MRIDoubleVec& v1, const MRIDoubleVec& v2, MRIDoubleVec& resVec){
  resVec[0] = v1[1] * v2[2] - v2[1] * v1[2];
  resVec[1] = v1[2] * v2[0] - v2[2] * v1[0];
  resVec[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

// ============================
// CHECK IF POINT IS INSIDE BOX
// ============================
bool isPointInsideBox(double xCoord, double yCoord, double zCoord, const MRIDoubleVec& limitBox){
    return ((xCoord>=limitBox[0])&&(xCoord<=limitBox[1]))&&
           ((yCoord>=limitBox[2])&&(yCoord<=limitBox[3]))&&
           ((zCoord>=limitBox[4])&&(zCoord<=limitBox[5]));
}

// =======================
// PRINT BIN ARRAY TO FILE
// =======================
void printBinArrayToFile(string fileName, int numberOfBins, const MRIDoubleVec& binCenter, const MRIDoubleVec& binValues){
  // Open Output File
	FILE* outFile;
	outFile = fopen(fileName.c_str(),"w");
	// Write Header
  for(int loopA=0;loopA<numberOfBins;loopA++){
    fprintf(outFile,"%e %e\n",binCenter[loopA],binValues[loopA]);
  }
	// Close Output file
	fclose(outFile);  
}

// ======================
// APPLY LIMIT BOX FACTOR
// ======================
void applyLimitBoxFactors(double xFactor, double yFactor, double zFactor, MRIDoubleVec& limitBox){
  double center[3] = {0.0};
  double domainSize[3] = {0.0};
  center[0] = 0.5 * (limitBox[0] + limitBox[1]);
  center[1] = 0.5 * (limitBox[2] + limitBox[3]);
  center[2] = 0.5 * (limitBox[4] + limitBox[5]);
  domainSize[0] = (limitBox[1] - limitBox[0]);
  domainSize[1] = (limitBox[3] - limitBox[2]);
  domainSize[2] = (limitBox[5] - limitBox[4]);
  limitBox[0] = center[0] - 0.5*xFactor*domainSize[0];
  limitBox[1] = center[0] + 0.5*xFactor*domainSize[0];
  limitBox[2] = center[1] - 0.5*yFactor*domainSize[1];
  limitBox[3] = center[1] + 0.5*yFactor*domainSize[1];
  limitBox[4] = center[2] - 0.5*zFactor*domainSize[2];
  limitBox[5] = center[2] + 0.5*zFactor*domainSize[2];
}

// ====================
// PRINT MATRIX TO FILE
// ====================
void printMatrixToFile(string fileName, const MRIDoubleMat& Mat){
  int totalRows = Mat.size();
  int totalCols = Mat[0].size();
    // Open Output File
	FILE* outFile;
	outFile = fopen(fileName.c_str(),"w");
	// Write Header
  for(int loopA=0;loopA<totalRows;loopA++){
    for(int loopB=0;loopB<totalCols;loopB++){
      if (loopB<(totalCols-1)){
        if (fabs(Mat[loopA][loopB])>kMathZero) fprintf(outFile,"%15.6e,",Mat[loopA][loopB]);
        else fprintf(outFile,"%15.6e,",0.0);
      }else{
        if (fabs(Mat[loopA][loopB])>kMathZero) fprintf(outFile,"%15.6e",Mat[loopA][loopB]);
        else fprintf(outFile,"%15.6e,",0.0);
      }
    }
    fprintf(outFile,"\n");
  }
	// Close Output file
	fclose(outFile);  
}

// ==================
// READ LIST OF FILES
// ==================
void readFileList(std::string listName, std::vector<std::string> &fileList){
  std::string buffer;
  // Assign File
  std::ifstream lFile;
  lFile.open(listName.c_str());
  // Loop through the File
  while (std::getline(lFile,buffer)){
    fileList.push_back(buffer);
  }
  // Close File
  lFile.close();
}

// =====================
// READ MATRIX FROM FILE
// =====================
void readMatrixFromFile(std::string inFileName,int& nrow,int& ncol,std::vector<std::vector<double> > &inMat){
  // Open File
  std::ifstream inFile;
  std::string buffer;
  std::vector<std::string> strs;
  inFile.open(inFileName.c_str());
  if (inFile.is_open()){
    while (!inFile.eof()){
      // Read One Line of Input File
      std::getline(inFile,buffer);
      if(buffer != ""){
        nrow++;
        // Get Number of Row
        if(nrow == 1){
          boost::split(strs,buffer,boost::is_any_of(" "));
          ncol = strs.size();
        }
      }
    }
    // ALLOCATE
    inMat.resize(nrow);
    for(int loopA = 0;loopA<nrow;loopA++){
      inMat[loopA].resize(ncol);
    }
    int count = -1;
    inFile.clear();
    inFile.seekg(0, std::ios::beg);
    while (!inFile.eof()){
      // Read One Line of Input File
      std::getline(inFile,buffer);
      if(buffer != ""){
        count++;
        // Get Number of Row
        boost::split(strs,buffer,boost::is_any_of(" "));
        for(unsigned int loopA=0;loopA<strs.size();loopA++){
          inMat[count][loopA] = atof(strs[loopA].c_str());
        }
      }
    }
  }
  inFile.close();
}

// ====================
// GET MEDIAN OF VECTOR
// ====================
double getMedian(MRIDoubleVec& v){
  size_t n = v.size() / 2;
  nth_element(v.begin(), v.begin()+n, v.end());
  return v[n];
}

// ========
// GET MEAN
// ========
double getMean(const MRIDoubleVec& v){
  double av = 0.0;
  for(size_t loopA=0;loopA<v.size();loopA++){
    av = av + v[loopA];
  }
  av = (av/((double)v.size()));
  return av;
}

// =============================
// GET EIGENVALUES OF 3x3 MATRIX
// =============================
void compute3x3MatrixEigenvals(double A[3][3], double root[3]){
  double inv3 = (1.0/3.0);
  double root3 = sqrt(3.0);
  double a00 = (double)A[0][0];
  double a01 = (double)A[0][1];
  double a02 = (double)A[0][2];
  double a11 = (double)A[1][1];
  double a12 = (double)A[1][2];
  double a22 = (double)A[2][2];
  double c0 = a00*a11*a22 + 2.0*a01*a02*a12 - a00*a12*a12 - a11*a02*a02 - a22*a01*a01;
  double c1 = a00*a11 - a01*a01 + a00*a22 - a02*a02 + a11*a22 - a12*a12;
  double c2 = a00 + a11 + a22;
  double c2Div3 = c2*inv3;
  double aDiv3 = (c1 - c2*c2Div3)*inv3;
  if (aDiv3 > 0.0){
    aDiv3 = 0.0;
  }
  double mbDiv2 = 0.5*(c0 + c2Div3*(2.0*c2Div3*c2Div3 - c1));
  double q = mbDiv2*mbDiv2 + aDiv3*aDiv3*aDiv3;
  if (q > 0.0) {
    q = 0.0;
  }
  double magnitude = sqrt(-aDiv3);
  double angle = atan2(sqrt(-q),mbDiv2)*inv3;
  double cs = cos(angle);
  double sn = sin(angle);
  root[0] = c2Div3 + 2.0*magnitude*cs;
  root[1] = c2Div3 - magnitude*(cs + root3*sn);
  root[2] = c2Div3 - magnitude*(cs - root3*sn);
  // SORT
  double temp = 0.0;
  if(root[0]<root[1]){
    temp = root[0];
    root[0] = root[1];
    root[1] = temp;
  }
  if(root[0]<root[2]){
    temp = root[0];
    root[0] = root[2];
    root[2] = temp;
  }
  if(root[1]<root[2]){
    temp = root[1];
    root[1] = root[2];
    root[2] = temp;
  }
}

// BOUBLE SORT ARRAY OF INTEGER
void sortIntArray(std::vector<int> &faceIds){
  int firstIdx = 0;
  int secondIdx = 0;
  int temp = 0;
  for(size_t loopA=0;loopA<faceIds.size()-1;loopA++){
    for(size_t loopB=loopA+1;loopB<faceIds.size();loopB++){
      firstIdx = faceIds[loopA];
      secondIdx = faceIds[loopB];
      if(firstIdx>secondIdx){
        temp = faceIds[loopA];
        faceIds[loopA] = faceIds[loopB];
        faceIds[loopB] = temp;
      }
    }
  }
}

// ===================================
// CHECK THAT TWO VECTORS ARE THE SAME
// ===================================
bool isSameIntVector(const MRIIntVec& one, const MRIIntVec& two){
  MRIIntVec first(one);
  MRIIntVec second(two);
  // SORT THE TWO VECTORS FIRST
  std::sort(first.begin(),first.end());
  std::sort(second.begin(),second.end());
  bool res = true;
  if((first.size() == 0)||(second.size() == 0)){
    return false;
  }
  if(first.size()!=second.size()){
    return false;
  }
  for(size_t loopA=0;loopA<first.size();loopA++){
    res = res && (first[loopA] == second[loopA]);
  }
  return res;
}

// =======================
// FIND HOW MANY INTERVALS
// =======================
int findHowMany(double distance, const MRIDoubleVec& lengths){
  bool found = false;
  int count = 0;
  double currDist = 0.0;
  while (!found){
    found = ((currDist + 1.0e-4) > distance);
    if(!found){
      currDist += lengths[count];
      count++;
    }
  }
  return count;
}

// ============================
// CHECK ERROR GENERATED BY MPI
// ============================
void checkMpiError(int mpiError){
  switch(mpiError){
    case MPI_SUCCESS:
      break;
    case MPI_ERR_COMM:
      printf("Invalid communicator. A common error is to use a null communicator in a call (not even allowed in MPI_Comm_rank).\n");
      exit(1);
      break;
    case MPI_ERR_COUNT:
      printf("Invalid count argument. Count arguments must be non-negative; a count of zero is often valid.\n");
      exit(1);
      break;
    case MPI_ERR_TYPE:
      printf("Invalid datatype argument. May be an uncommitted MPI_Datatype (see MPI_Type_commit).\n");
      exit(1);
      break;
    case MPI_ERR_BUFFER:
      printf("Invalid buffer pointer. Usually a null buffer where one is not valid.\n");
      exit(1);
      break;
    case MPI_ERR_ROOT:
      printf("Invalid root. The root must be specified as a rank in the communicator. Ranks must be between zero and the size of the communicator minus one.\n");
      exit(1);
      break;
  }
}

// =============================
// GET THRESHOLD QUANTITY STRING
// =============================
string getThresholdQtyString(int thresholdQty){
  string res;
  switch(thresholdQty){
    case kQtyConcentration:
      res = "Concentration";
      break;
    case kQtyVelocityX:
      res = "VelocityX";
      break;
    case kQtyVelocityY:
      res = "VelocityY";
      break;
    case kQtyVelocityZ:
      res = "VelocityZ";
      break;
  }
  return res;
}

// =========================
// GET THRESHOLD TYPE STRING
// =========================
string getThresholdTypeString(int thresholdType){
  string res;
  switch(thresholdType){
    case kCriterionLessThen:
      res = "Set to zero if less than threshold";
      break;
    case kCriterionGreaterThen:
      res = "Set to zero if greater than threshold";
      break;
    case kCriterionABSLessThen:
      res = "Set to zero if ABS less than threshold";
      break;
    case kCriterionABSGreaterThen:
      res = "Set to zero if ABS greater than threshold";
      break;
  }
  return res;
}

// =========================
// GET MIN OF STD INT VECTOR
// =========================
int getMinInt(vector<int> vec){
  int minVal = INT_MAX;
  for(size_t loopA=0;loopA<vec.size();loopA++){
    if(vec[loopA] < minVal){
      minVal = vec[loopA];
    }
  }
  return minVal;
}

// ===================================
// CHECK IF STRING IS A FLOATING POINT
// ===================================
bool isFloat(string token){
  std::istringstream iss(token);
  float f;
  iss >> noskipws >> f; // noskipws considers leading whitespace invalid
  // Check the entire string was consumed and if either failbit or badbit is set
  return iss.eof() && !iss.fail(); 
}

// ====================
// READ TABLE FROM FILE
// ====================
int readTableFromFile(string fileName, MRIDoubleMat& table, bool skipFirstRow){
  // Open File
  ifstream myReadFile;
  string buffer;
  int lineCount = 0;
  std::vector<std::string> tokens;
  std::vector<double> currParams;
  myReadFile.open(fileName.c_str());
  if (myReadFile.is_open()) {
    if(skipFirstRow){
      std::getline(myReadFile,buffer);
    }
    // Loop through the File
    while (!myReadFile.eof()){
      // Read One Line of Input File
      std::getline(myReadFile,buffer);
      if(!buffer.empty()){
        lineCount++;
        boost::split(tokens, buffer, boost::is_any_of(","));
        currParams.clear();
        for(int loopA=0;loopA<tokens.size();loopA++){
          try{
            currParams.push_back(atof(tokens[loopA].c_str()));
          }catch(...){
            printf("Invalid Table File. Exiting.\n");
            return 1;
          }
        }
        table.push_back(currParams);
      }
    }
  }else{
    printf("Cannot Open File. Exiting.\n");
    return 1;
  }
  // Close File
  myReadFile.close();
  // Return
  return 0;
}

// ===================
// PRINT DOUBLE MATRIX
// ===================
void printDoubleMatToFile(string fileName, const MRIDoubleMat& faceNormals){
  // Open Output File
  FILE* f;
  f = fopen(fileName.c_str(),"w");
  // Write Header
  for(int loopA=0;loopA<faceNormals.size();loopA++){
    for(int loopB=0;loopB<faceNormals[loopA].size();loopB++){
      fprintf(f,"%f ",faceNormals[loopA][loopB]);
    }
    fprintf(f,"\n");
  }
  // Close Output file
  fclose(f);
}

// ================
// PRINT INT MATRIX
// ================
void printIntMatToFile(string fileName, const MRIIntMat& mat){
  // Open Output File
  FILE* f;
  f = fopen(fileName.c_str(),"w");
  // Write Header
  for(int loopA=0;loopA<mat.size();loopA++){
    for(int loopB=0;loopB<mat[loopA].size();loopB++){
      fprintf(f,"%d ",mat[loopA][loopB]);
    }
    fprintf(f,"\n");
  }
  // Close Output file
  fclose(f);
}

// ===================
// PRINT DOUBLE VECTOR
// ===================
void printDoubleVecToFile(string fileName, const MRIDoubleVec& vec){
  // Open Output File
  FILE* f;
  f = fopen(fileName.c_str(),"w");
  // Write Header
  for(int loopA=0;loopA<vec.size();loopA++){
    fprintf(f,"%f\n",vec[loopA]);
  }
  // Close Output file
  fclose(f);
}

// ===================
// PRINT DOUBLE VECTOR
// ===================
void printIntVecToFile(string fileName, const MRIIntVec& vec){
  // Open Output File
  FILE* f;
  f = fopen(fileName.c_str(),"w");
  // Write Header
  for(int loopA=0;loopA<vec.size();loopA++){
    fprintf(f,"%d\n",vec[loopA]);
  }
  // Close Output file
  fclose(f);
}

// ===================
// PRINT DOUBLE MATRIX
// ===================
void printDoubleArrayToFile(string fileName, int size, double* vec){
  // Open Output File
  FILE* f;
  f = fopen(fileName.c_str(),"w");
  // Write Header
  for(int loopA=0;loopA<size;loopA++){
    fprintf(f,"%f\n",vec[loopA]);
  }
  // Close Output file
  fclose(f);
}

// ==============
// NORMALIZE BINS
// ==============
void normalizeBinArray(MRIDoubleVec& binArray,double currInterval){
  double sum = 0.0;
  // Compute Summation
  for(int loopA=0;loopA<binArray.size();loopA++){
    sum += binArray[loopA];
  }
  // Exit if zero sum
  if (fabs(sum)<kMathZero){
    return;
  }
  // Normalize
  for(int loopA=0;loopA<binArray.size();loopA++){
    binArray[loopA] /= (sum*currInterval);
  }  
}

// =============
// ASSIGN TO BIN
// =============
void assignToBin(double currValue, int numberOfBins, const MRIDoubleVec& binMin, const MRIDoubleVec& binMax, MRIDoubleVec& binArray){
  bool found = false;
  int count = 0;
  bool isMoreThanMin = false;
  bool isLessThanMax = false;
  while ((!found)&&(count<numberOfBins)){
    if (fabs(currValue-binMin[0])<kMathZero){
      isMoreThanMin = (currValue >= binMin[count] - kMathZero);
    }else{
      isMoreThanMin = (currValue > binMin[count]);
    }
    if (fabs(currValue-binMin[numberOfBins-1])<kMathZero){
      isLessThanMax = (currValue <= binMax[count]);
    }else{
      isLessThanMax = (currValue <= binMax[count] + kMathZero);
    }
    found = (isMoreThanMin)&&(isLessThanMax);
    // Update
    if (!found){
      count++;
    }
  }
  if (found){
    // Increase Bin Count
    binArray[count] = binArray[count] + 1.0;
  }else{
    throw MRIException("Error: Value Cannot fit in Bin.\n");
  } 
}

}
