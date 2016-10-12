#ifndef MRIUTILS_H
#define MRIUTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
// String Utilities
#include <boost/algorithm/string.hpp>
// MPI
#include "mpi.h"

#include "mriConstants.h"
#include "mriCoordItem.h"
#include "mriStreamline.h"
#include "mriThresholdCriteria.h"
#include "mriTypes.h"
#include "schMessages.h"

using namespace std;

// MRI UTILITIES
namespace MRIUtils{

  // GLOBAL VARIABLE FOR BOOST RANDOM NUMBER GENERATION
  static boost::random::mt19937 realGen;

  // Write Help
  void writeProgramHelp();

  // Home made Int to String Conversion
  string intToStr(int number);

  // Home made Float to String Conversion
  string floatToStr(double number);

  // Own String Split Method
  MRIStringVec extractSubStringFromBufferMS(string Buffer);

  // Compute Euclidean Norm
  template<typename Type> Type do3DEucNorm(vector<Type> v);

  // Normalize 3-component vector
  void normalize3DVector(MRIDoubleVec& v);

  // Insert Value in String 
  template<typename Type> void insertInList(Type value, vector<Type>& coords);

  // Sample from Standard Gaussian Distribution
  double generateStandardGaussian(double stDev);

  // Round Floating Point Value
  int round(double value);

  // Write Abcissae and Ordinates to file
  void writeGraphToFile(string fileName, int vecSize, MRIDoubleVec& vecX, MRIDoubleVec& vecY);

  // Perform 3D External Product
  void do3DExternalProduct(const MRIDoubleVec& v1, const MRIDoubleVec& v2, MRIDoubleVec& resVec);

  // Check if Point is Inside/Outside Box
  bool isPointInsideBox(double xCoord, double yCoord, double zCoord, const MRIDoubleVec& limitBox);

  // Print Binary Array To File
  void printBinArrayToFile(string fileName, int numberOfBins, const MRIDoubleVec& binCenter, const MRIDoubleVec& binValues);

  // Apply Limit Box Factor
  void applyLimitBoxFactors(double xFactor, double yFactor, double zFactor, MRIDoubleVec& limitBox);

  // Print Matrix To File
  void printMatrixToFile(string fileName, int totalRows, int totalCols, MRIDoubleMat& Mat);

  // Write Streamlines to VTK
  void printStreamlinesToVTK(vector<MRIStreamline*> streamlines,std::string fileName);

  // Read Streamlines From VTK
  void readStreamlinesFromLegacyVTK(string fileName, std::vector<MRIStreamline*> &streamlines);

  // Read List of Files
  void readFileList(string listName, MRIStringVec& fileList);

  // Read Matrix From File
  void readMatrixFromFile(string inFileName, int& nrow, int& ncol, MRIDoubleMat& inMat);

  // Get Median of Double Vector
  double GetMedian(const MRIDoubleVec& v);

  // Get Mean of Double Vector
  double GetMean(const MRIDoubleVec& v);

  // Get Eigenvalues of 3x3 Matrix
  void compute3x3MatrixEigenvals(double A[3][3], double root[3]);

  // Bouble Sort Array of Integer
  void sortIntArray(MRIIntVec& faceIds);

  // Check that two integer vectors are the same
  bool isSameIntVector(const MRIIntVec& first, const MRIIntVec& second);

  // FIND HOW MANY INTERVALS
  int findHowMany(double distance, std::vector<double> lengths);

  // CHECK ERROR GENERATED BY MPI
  void checkMpiError(int mpiError);

  // GET THRESHOLD QUANTITY STRING
  string getThresholdQtyString(int thresholdQty);

  // GET THRESHOLD TYPE STRING
  string getThresholdTypeString(int thresholdType);

  // GET MIN OF STD INT VECTOR
  int getMinInt(vector<int> vec);

  // CHECK IF STRING IS A FLOATING POINT
  bool isFloat(string token);

  // READ TABLE FROM FILE
  int readTableFromFile(string fileName, MRIDoubleMat& table,bool skipFirstRow);

}
#endif //MRIUTILS_H
