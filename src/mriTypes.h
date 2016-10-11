#ifndef MRITYPES_H
#define MRITYPES_H

# include <vector>

// ARRAYS
typedef std::vector< std::vector<double> > MRIDoubleMat;
typedef std::vector<double> MRIDoubleVec;
typedef std::vector< std::vector<int> > MRIIntMat;
typedef std::vector<int> MRIIntVec;
typedef std::vector<bool> MRIBoolVec;
typedef std::vector<std::string> MRIStringVec;

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

