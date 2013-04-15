#ifndef MRIUTILS_H
#define MRIUTILS_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sstream>
#include <vector>
// Random Number Generator
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

#include "mriConstants.h"
#include "mriCoordItem.h"
#include "mriStreamline.h"

// MRI UTILITIES
namespace MRIUtils{
  // GLOBAL VARIABLE FOR GENERATION
  //static boost::random::mt19937 intGen;
  //static boost::random::mt19937 realGen;
  //static boost::random::mt19937 normGen;
  
// Write Help
inline void WriteProgramHelp()
{
   printf("MRIApp v0.0.4 - Daniele Schiavazzi - 2012\n");
   printf("usage MRIApp [options] [files]||[sequence]\n");
};	
// Compare 4 Integers
inline bool Compare4Integer4(int Int1,int Int2,int Int3,int Int4)
{
  bool myresult = true;
  if ((Int1!=Int2)||(Int1!=Int3)||(Int1!=Int4)||(Int2!=Int3)||(Int2!=Int4)||(Int3!=Int4))
  { 
	myresult = false; 
  }
  return myresult;
};
// Compare 4 float
inline bool Compare4Single(float Sing1, float Sing2, float Sing3, float Sing4)
{
  bool result = true;
  if ((fabs(Sing1-Sing2)>kMathZero)||
     (fabs(Sing1-Sing3)>kMathZero)||
     (fabs(Sing1-Sing4)>kMathZero)||
     (fabs(Sing2-Sing3)>kMathZero)||
     (fabs(Sing2-Sing4)>kMathZero)||
     (fabs(Sing3-Sing4)>kMathZero)) result = false;
	return result;
};
// Convert Integer to String
inline std::string IntToStr(int number)
{
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
};
// Convert Double to String
inline std::string FloatToStr(double number)
{
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
};

inline std::vector<std::string> ExctractSubStringFromBufferMS(std::string Buffer)
{
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
};

// Eval 2 Norm of vector
template<typename Type>
inline Type Do3DEucNorm(Type* v)
{
  Type norm2 = 0.0;
  for(int i = 0;i<kNumberOfDimensions;i++) {
    norm2 += v[i]*v[i];
  }  
  return sqrt(norm2);
};

// ===================
// NORMALIZE 3D VECTOR
// ===================
inline void Normalize3DVector(double* v)
{
  double norm = Do3DEucNorm(v);
  if(norm>kMathZero){
    for(int LoopA=0;LoopA<kNumberOfDimensions;LoopA++)
    {
      v[LoopA] = v[LoopA]/norm;
    }
  }
};

// =======================================
// INSERTION IN VECTORS WITH NO DUPLICATES
// =======================================

// INSERT IN DOUBLE LIST
inline void InsertInIntList(int Item,int &TotalFaces,std::vector<int> &FacesID)
{
  for(int LoopA=0;LoopA<TotalFaces;LoopA++)	
	{
		if (FacesID[LoopA]==Item) return;
	}	
  TotalFaces++;
	FacesID.push_back(Item);
};

// INSERT IN INT LIST
inline void InsertInDoubleList(double Value,int &TotalCoords,std::vector<double> &Coords)
{
  double distance;	  
  for(int LoopA=0;LoopA<TotalCoords;LoopA++)
  { 
    distance = fabs(Coords[LoopA]-Value);
    if (distance<kMathZero) return;
  }
  // Insert New Value
  TotalCoords++;
  Coords.push_back(Value);
};

// Generate Uniform Integer (global Variable Defined)
inline int GenerateUniformIntegers(int lowIdx, int upIdx) {
    boost::random::uniform_int_distribution<> dist(lowIdx, upIdx);
    //return dist(intGen);
    return 0;
}

// GENERATE RANDOM VECTOR
inline double* GenerateUniform01RandomVector(double& generator){
  double* resVec = new double[3];
  boost::random::uniform_real_distribution<> dist(0.0,1.0);
  // Generate Numbers Uniformly
  for(int loopA=0;loopA<3;loopA++){
		//resVec[loopA] = dist(realGen);
	}
  // Normalize
  MRIUtils::Normalize3DVector(resVec);
	return resVec;
}

// PERTURB ITS COORDINATES WITH GAUSSIAN NOISE
inline double GenerateStandardGaussian(double stDev){
  // Allocate Vector
  boost::random::normal_distribution<> dist(0.0,stDev);
  // Add Random Component
  //return dist(realGen);
  return 0.0;
}


// PERTURB ITS COORDINATES WITH GAUSSIAN NOISE
inline void PerturbVectorGaussian(int vecSize, MRICoordItem* coordItemArray, double stDev){
  // Allocate Vector
  boost::random::normal_distribution<> dist(0.0,stDev);
  // Add Random Component
  for(int loopA=0;loopA<vecSize;loopA++){
    //coordItemArray[loopA].x += dist(normGen);
    //coordItemArray[loopA].y += dist(normGen);
    //coordItemArray[loopA].z += dist(normGen);
  }
}

// Round Values
inline int round(double value){
  return (int)floor(value+0.5);
}

// Write Vector Graph To File
inline void WriteGraphToFile(std::string fileName, int vecSize, std::vector<double> &vecX, std::vector<double> &vecY){
  // Open Output File
	FILE* outFile;
	outFile = fopen(fileName.c_str(),"w");
	// Write Header
  for(int loopA=0;loopA<vecSize;loopA++){
    fprintf(outFile,"%d %e %e\n",loopA,vecX[loopA],vecY[loopA]);
  }
	// Close Output file
	fclose(outFile);  
}

// PERFORM EXTERNAL PRODUCT
inline void Do3DExternalProduct(double* v1, double* v2, double* resVec){
  resVec[0] = v1[1] * v2[2] - v2[1] * v1[2];
  resVec[1] = v1[2] * v2[0] - v2[2] * v1[0];
  resVec[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

// CHECK IF POINT IS INSIDE BOX
inline bool IsPointInsideBox(double xCoord, double yCoord, double zCoord, double* limitBox){
    return ((xCoord>=limitBox[0])&&(xCoord<=limitBox[1]))&&
           ((yCoord>=limitBox[2])&&(yCoord<=limitBox[3]))&&
           ((zCoord>=limitBox[4])&&(zCoord<=limitBox[5]));
}

// PRINT BIN ARRAY TO FILE
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

// APPLY LIMIT BOX FACTOR
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

// PRINT MATRIX TO FILE
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

// ================
// Print Streamlines to VTK
// ================
inline void PrintStreamlinesToVTK(std::vector<MRIStreamline*> streamlines,std::string fileName){
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
  fprintf(outFile,"CELLS %d %d\n",streamlines.size(),streamlines.size()+totalSLPoints);

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
  fprintf(outFile,"CELL_TYPES %d\n",streamlines.size());
  for (unsigned int loopA=0;loopA<streamlines.size();loopA++){
    fprintf(outFile,"4\n");
  }

}

#endif //MRIUTILS_H
