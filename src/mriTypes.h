#ifndef MRITYPES_H
#define MRITYPES_H

# include <vector>
# include <string>

using namespace std;

// ARRAYS
typedef vector< vector<double> > MRIDoubleMat;
typedef vector<double> MRIDoubleVec;
typedef vector< vector<int> > MRIIntMat;
typedef vector<int> MRIIntVec;
typedef vector<bool> MRIBoolVec;
typedef vector< string > MRIStringVec;

// RELATIVE POSITION BETWEEN EDGE AND FACE
enum EdgeFacePositionType{ptTop,ptBottom,ptLeft,ptRight};

// AUXILIARY FACE
struct mriFace{
  int number;
  std::vector<int> connections;
};

// AUXILIARY EDGE
struct mriEdge{
  int number;
  std::vector<int> connections;
};

// ===================
// TYPES FOR PLT FILES
// ===================
enum pltFileTypes{
  pltUNIFORM,
  pltSTRUCTURED
};

struct PLTOptionRecord{
  int i;
  int j;
  int k;
  int N;
  int E;
  pltFileTypes type;
};

// ===================
// TYPES FOR PLT FILES
// ===================
struct vtkStructuredPointsOptionRecord{
  bool isASCII;
  bool isValidDataset;
  int dimensions[3];
  double origin[3];
  double spacing[3];
  int numDefined;
  bool isDefined[5];
  MRIIntVec dataBlockStart;
  MRIIntVec dataBlockType;
  MRIBoolVec dataBlockRead;
};

#endif // MRITYPES_H

