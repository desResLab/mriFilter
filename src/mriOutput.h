#ifndef MRIOUTPUT_H
#define MRIOUTPUT_H

# include <vector>
# include <string>

# include "mriTypes.h"

using namespace std;

class MRIOutput{
  public:
    // DATA MEMBERS
    string name;
    int totComponents;
    MRIDoubleVec values;

    // CONSTRUCTOR
    MRIOutput(string locName, int comp);
};

#endif // MRIOUTPUT_H
