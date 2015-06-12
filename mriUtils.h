#ifndef MRIUTILS_H
#define MRIUTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
// String Utilities
#include <boost/algorithm/string.hpp>
#include <regex>
// MPI
#include "mpi.h"

#include "mriConstants.h"
#include "mriCoordItem.h"
#include "mriStreamline.h"
#include "mriThresholdCriteria.h"
#include "schMessages.h"

using namespace std;

// MRI UTILITIES
namespace MRIUtils{
  // GLOBAL VARIABLE FOR GENERATION
  //static boost::random::mt19937 intGen;
  static boost::random::mt19937 realGen;
  //static boost::random::mt19937 normGen;

// ==================
// WRITE PROGRAM HELP
// ==================
inline void WriteProgramHelp(){
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

// ==================
// COMPARE 4 INTEGERS
// ==================
inline bool Compare4Integer4(int Int1,int Int2,int Int3,int Int4){
  bool myresult = true;
  if ((Int1!=Int2)||(Int1!=Int3)||(Int1!=Int4)||(Int2!=Int3)||(Int2!=Int4)||(Int3!=Int4))
  { 
	myresult = false; 
  }
  return myresult;
}

// ================
// COMPARE 4 FLOATS
// ================
inline bool Compare4Single(float Sing1, float Sing2, float Sing3, float Sing4){
  bool result = true;
  if ((fabs(Sing1-Sing2)>kMathZero)||
     (fabs(Sing1-Sing3)>kMathZero)||
     (fabs(Sing1-Sing4)>kMathZero)||
     (fabs(Sing2-Sing3)>kMathZero)||
     (fabs(Sing2-Sing4)>kMathZero)||
     (fabs(Sing3-Sing4)>kMathZero)) result = false;
	return result;
}

// =========================
// CONVERT INTEGER TO STRING
// =========================
inline std::string IntToStr(int number){
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

// ========================
// CONVERT DOUBLE TO STRING
// ========================
inline std::string FloatToStr(double number){
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

// ===============================================
// EXTRACT STRING FROM BUFFER WITH MULTIPLE SPACES
// ===============================================
inline std::vector<std::string> ExctractSubStringFromBufferMS(std::string Buffer){
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
template<typename Type>
inline Type Do3DEucNorm(Type* v){
  Type norm2 = 0.0;
  for(int i = 0;i<kNumberOfDimensions;i++){
    norm2 += v[i]*v[i];
  }  
  return sqrt(norm2);
}
template<typename Type>
inline Type Do3DEucNorm(std::vector<Type> v){
  Type norm2 = 0.0;
  for(int i = 0;i<kNumberOfDimensions;i++){
    norm2 += v[i]*v[i];
  }
  return sqrt(norm2);
}


// ===================
// NORMALIZE 3D VECTOR
// ===================
inline void Normalize3DVector(double* v){
  double norm = Do3DEucNorm(v);
  if(norm>kMathZero){
    for(int LoopA=0;LoopA<kNumberOfDimensions;LoopA++)
    {
      v[LoopA] = v[LoopA]/norm;
    }
  }
}
inline void Normalize3DVector(std::vector<double> &v){
  double norm = Do3DEucNorm(v);
  if(norm>kMathZero){
    for(int LoopA=0;LoopA<kNumberOfDimensions;LoopA++)
    {
      v[LoopA] = v[LoopA]/norm;
    }
  }
}


// =======================================
// INSERTION IN VECTORS WITH NO DUPLICATES
// =======================================

// INSERT IN DOUBLE LIST
inline void InsertInIntList(int Item,int &TotalFaces,std::vector<int> &FacesID){
  for(int LoopA=0;LoopA<TotalFaces;LoopA++){
    if (FacesID[LoopA]==Item) return;
  }
  TotalFaces++;
  FacesID.push_back(Item);
}

// INSERT IN INT LIST
inline void InsertInDoubleList(double Value,int &TotalCoords,std::vector<double> &Coords){
  double distance;	  
  for(int LoopA=0;LoopA<TotalCoords;LoopA++)
  { 
    distance = fabs(Coords[LoopA]-Value);
    if (distance<kMathZero) return;
  }
  // Insert New Value
  TotalCoords++;
  Coords.push_back(Value);
}

// ==================================================
// Generate Uniform Integer (global Variable Defined)
// ==================================================
/*inline int GenerateUniformIntegers(int lowIdx, int upIdx) {
    boost::random::uniform_int_distribution<> dist(lowIdx, upIdx);
    //return dist(intGen);
    return 0;
}*/

// ======================
// GENERATE RANDOM VECTOR
// ======================
/*inline double* GenerateUniform01RandomVector(double& generator){
  double* resVec = new double[3];
  boost::random::uniform_real_distribution<> dist(0.0,1.0);
  // Generate Numbers Uniformly
  for(int loopA=0;loopA<3;loopA++){
		//resVec[loopA] = dist(realGen);
	}
  // Normalize
  MRIUtils::Normalize3DVector(resVec);
	return resVec;
}*/

// ==========================
// GENERATE STANDARD GAUSSIAN
// ==========================
inline double GenerateStandardGaussian(double stDev){
  // Allocate Vector
  boost::random::normal_distribution<> dist(0.0,stDev);
  // Add Random Component
  return dist(realGen);
}

// ===========================================
// PERTURB ITS COORDINATES WITH GAUSSIAN NOISE
// ===========================================
/*inline void PerturbVectorGaussian(int vecSize, MRICoordItem* coordItemArray, double stDev){
  // Allocate Vector
  boost::random::normal_distribution<> dist(0.0,stDev);
  // Add Random Component
  for(int loopA=0;loopA<vecSize;loopA++){
    //coordItemArray[loopA].x += dist(normGen);
    //coordItemArray[loopA].y += dist(normGen);
    //coordItemArray[loopA].z += dist(normGen);
  }
}*/

// ============
// ROUND VALUES
// ============
inline int round(double value){
  return (int)floor(value+0.5);
}

// ==========================
// WRITE VECTOR GRAPH TO FILE
// ==========================
inline void WriteGraphToFile(std::string fileName, int vecSize, std::vector<double> &vecX, std::vector<double> &vecY){
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
inline void Do3DExternalProduct(double* v1, double* v2, double* resVec){
  resVec[0] = v1[1] * v2[2] - v2[1] * v1[2];
  resVec[1] = v1[2] * v2[0] - v2[2] * v1[0];
  resVec[2] = v1[0] * v2[1] - v2[0] * v1[1];
}
inline void Do3DExternalProduct(std::vector<double> v1, std::vector<double> v2, std::vector<double> &resVec){
  resVec[0] = v1[1] * v2[2] - v2[1] * v1[2];
  resVec[1] = v1[2] * v2[0] - v2[2] * v1[0];
  resVec[2] = v1[0] * v2[1] - v2[0] * v1[1];
}


// ============================
// CHECK IF POINT IS INSIDE BOX
// ============================
inline bool IsPointInsideBox(double xCoord, double yCoord, double zCoord, double* limitBox){
    return ((xCoord>=limitBox[0])&&(xCoord<=limitBox[1]))&&
           ((yCoord>=limitBox[2])&&(yCoord<=limitBox[3]))&&
           ((zCoord>=limitBox[4])&&(zCoord<=limitBox[5]));
}

// =======================
// PRINT BIN ARRAY TO FILE
// =======================
inline void PrintBinArrayToFile(std::string fileName, int numberOfBins, double* binCenter, double* binValues){
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
inline void ApplylimitBoxFactors(double xFactor, double yFactor, double zFactor,double* limitBox){
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
inline void PrintMatrixToFile(std::string fileName, int totalRows, int totalCols, double** Mat){
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

// ========================
// PRINT STREAMLINES TO VTK
// ========================
inline void PrintStreamlinesToVTK(std::vector<MRIStreamline*> streamlines,std::string fileName){
  // Write Message
  WriteSchMessage(std::string("\n"));
  WriteSchMessage(std::string("Writing Streamlines to VTK file...\n"));
  // Open Output File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  // Write Header
  fprintf(outFile,"# vtk DataFile Version 2.0\n");
  fprintf(outFile,"Streamline Data\n");
  fprintf(outFile,"ASCII\n");

  // Writre Data Set
  fprintf(outFile,"DATASET UNSTRUCTURED_GRID\n");

  // Eval the total Number Of Points
  int totalSLPoints = 0;
  for (unsigned int loopA=0;loopA<streamlines.size();loopA++){
    totalSLPoints += streamlines[loopA]->totalPoints;
  }
  // Write Point Header
  fprintf(outFile,"POINTS %d float\n",totalSLPoints);
  // Write Point Coordinates
  for (unsigned int loopA=0;loopA<streamlines.size();loopA++){
    for (int loopB=0;loopB<streamlines[loopA]->totalPoints;loopB++){
      fprintf(outFile,"%e %e %e\n",streamlines[loopA]->xCoords[loopB],streamlines[loopA]->yCoords[loopB],streamlines[loopA]->zCoords[loopB]);
    }
  }

  // Write Cells Header
  fprintf(outFile,"CELLS %d %d\n",(int)streamlines.size(),(int)(streamlines.size()+totalSLPoints));

  // Write Cells Definition
  int nodeCount =0;
  for (unsigned int loopA=0;loopA<streamlines.size();loopA++){
    fprintf(outFile,"%d ",streamlines[loopA]->totalPoints);
    for (int loopB=0;loopB<streamlines[loopA]->totalPoints;loopB++){
      fprintf(outFile,"%d ",nodeCount);
      nodeCount++;
    }
    fprintf(outFile,"\n");
  }

  // Write Cells Types Header
  fprintf(outFile,"CELL_TYPES %d\n",(int)streamlines.size());
  for (unsigned int loopA=0;loopA<streamlines.size();loopA++){
    fprintf(outFile,"4\n");
  }

  // Close Output file
  fclose(outFile);
}

// =========================
// READ STREAMLINES FROM VTK
// =========================
inline void ReadStreamlinesFromLegacyVTK(std::string fileName, std::vector<MRIStreamline*> &streamlines){
  // Write Message
  WriteSchMessage(std::string("\n"));
  WriteSchMessage(std::string("Reading Streamlines from file...\n"));
  // Declare input File
  std::ifstream infile;
  infile.open(fileName.c_str());

  // Read Data From File
  std::string buffer;
  std::vector<std::string> tokenizedString;
  std::vector<double> xCoords;
  std::vector<double> yCoords;
  std::vector<double> zCoords;
  std::vector<int> idx;

  // Read Points
  bool doReadPoints = false;
  bool doReadLines = false;
  int totPoints = 0;
  int totLines = 0;
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
    // Tokenize String
    boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);

    // Parse tokenized String
    if (boost::to_upper_copy(tokenizedString[0]) == std::string("POINTS")){
      // Get total number of points
      totPoints = atoi(tokenizedString[1].c_str());
      // Activate flags
      doReadPoints = true;
    }else if (boost::to_upper_copy(tokenizedString[0]) == std::string("LINES")){
      // Get total number of Lines
      totLines = atoi(tokenizedString[1].c_str());
      // Activate flags
      doReadLines = true;
    }

    // Read Nodes
    int numPointsReadInLine = 0;
    if(doReadPoints){
      // Write Message
      WriteSchMessage(std::string("Reading Streamlines Points...\n"));
      int readPoints = 0;
      while(readPoints<totPoints){
        // Read line in input File
        std::getline(infile,buffer);
        // Trim String
        boost::trim(buffer);
        // Tokenize String
        boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
        // Count number of points read
        numPointsReadInLine = (int)floor(tokenizedString.size()/(double)3.0);
        // Increment Counter
        readPoints += numPointsReadInLine;
        // Add Points to list
        for(int loopA=0;loopA<numPointsReadInLine;loopA++){
          xCoords.push_back(atof(tokenizedString[loopA*3].c_str()));
          yCoords.push_back(atof(tokenizedString[loopA*3+1].c_str()));
          zCoords.push_back(atof(tokenizedString[loopA*3+2].c_str()));
        }
      }
      // Reset flag
      doReadPoints = false;
    }

    // Read Lines
    if(doReadLines){
      // Write Message
      WriteSchMessage(std::string("Reading Streamlines Lines...\n"));
      int readLines = 0;
      while(readLines<totLines){
        // Clear Container
        idx.clear();
        // Read line in input File
        std::getline(infile,buffer);
        // Trim String
        boost::trim(buffer);
        // Tokenize String
        boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
        // Create New Streamline
        MRIStreamline* sl = new MRIStreamline();
        // Add Points to list
        for(unsigned int loopA=1;loopA<tokenizedString.size();loopA++){
          sl->totalPoints++;
          sl->xCoords.push_back(xCoords[atoi(tokenizedString[loopA].c_str())]);
          sl->yCoords.push_back(yCoords[atoi(tokenizedString[loopA].c_str())]);
          sl->zCoords.push_back(zCoords[atoi(tokenizedString[loopA].c_str())]);
        }
        streamlines.push_back(sl);
        // Update Lines Count
        readLines++;
      }
      // Reset flag
      doReadLines = false;
    }
  }
  // Close File
  infile.close();
}

// ==================
// READ LIST OF FILES
// ==================
inline void ReadFileList(std::string listName, std::vector<std::string> &fileList){
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
inline void ReadMatrixFromFile(std::string inFileName,int& nrow,int& ncol,std::vector<std::vector<double> > &inMat){
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
inline double GetMedian(std::vector<double> &v){
  size_t n = v.size() / 2;
  nth_element(v.begin(), v.begin()+n, v.end());
  return v[n];
}

// ========
// GET MEAN
// ========
inline double GetMean(std::vector<double> &v){
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
inline void Compute3x3MatrixEigenvals(double A[3][3], double root[3]){
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
inline void SortIntArray(std::vector<int> &faceIds){
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
inline bool isSameIntVector(std::vector<int> first, std::vector<int> second){
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
inline int FindHowMany(double distance, std::vector<double> lengths){
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

// =========================
// FIND IF A STRING IS FLOAT
// =========================
inline int isFloat(std::string currString){
  std::regex e("\\s*[+-]?([0-9]+\\.[0-9]*([Ee][+-]?[0-9]+)?|\\.[0-9]+([Ee][+-]?[0-9]+)?|[0-9]+[Ee][+-]?[0-9]+)");
  return std::regex_match(currString,e);
}

// ============================
// CHECK ERROR GENERATED BY MPI
// ============================
inline void checkMpiError(int mpiError){
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
inline string getThresholdQtyString(int thresholdQty){
  switch(thresholdQty){
    case kQtyPositionX:
      return string("PositionX");
      break;
    case kQtyPositionY:
      return string("PositionY");
      break;
    case kQtyPositionZ:
      return string("PositionZ");
      break;
    case kQtyConcentration:
      return string("Concentration");
      break;
    case kQtyVelocityX:
      return string("VelocityX");
      break;
    case kQtyVelocityY:
      return string("VelocityY");
      break;
    case kQtyVelocityZ:
      return string("VelocityZ");
      break;
  }
}

// =========================
// GET THRESHOLD TYPE STRING
// =========================
inline string getThresholdTypeString(int thresholdType){
  switch(thresholdType){
    case kCriterionLessThen:
      return string("Set to zero if less than threshold");
      break;
    case kCriterionGreaterThen:
      return string("Set to zero if greater than threshold");
      break;
    case kCriterionABSLessThen:
      return string("Set to zero if ABS less than threshold");
      break;
    case kCriterionABSGreaterThen:
      return string("Set to zero if ABS greater than threshold");
      break;
  }
}

}
#endif //MRIUTILS_H
