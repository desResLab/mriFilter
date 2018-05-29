#ifndef MRITYPES_H
#define MRITYPES_H

# include <vector>
# include <string>

using namespace std;

// ARRAYS
typedef vector< vector<double> > mriDoubleMat;
typedef vector<double>           mriDoubleVec;
typedef vector< vector<int> >    mriIntMat;
typedef vector<int>              mriIntVec;
typedef vector<bool>             mriBoolVec;
typedef vector< string >         mriStringVec;

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

struct pltOptionRecord{
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
  mriIntVec dataBlockStart;
  mriIntVec dataBlockType;
  mriBoolVec dataBlockRead;
};

#endif // MRITYPES_H

