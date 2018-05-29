#ifndef MRIOUTPUT_H
#define MRIOUTPUT_H

# include <vector>
# include <string>

# include "mriTypes.h"

using namespace std;

class mriOutput{
  public:
    // DATA MEMBERS
    string name;
    int totComponents;
    mriDoubleVec values;

    // CONSTRUCTOR
    mriOutput(string locName, int comp);
};

#endif // MRIOUTPUT_H
