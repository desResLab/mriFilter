#ifndef MRIOUTPUT_H
#define MRIOUTPUT_H

# include <vector>
# include <string>

using namespace std;

class MRIOutput
{
public:
    // DATA MEMBERS
    string name;
    int totComponents;
    vector<double> values;

    // CONSTRUCTOR
    MRIOutput(string locName, int comp);
};

#endif // MRIOUTPUT_H
