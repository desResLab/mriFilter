#ifndef MRIVOLDATA_H
#define MRIVOLDATA_H

class MRIVolData
{
  public:
    // Data Members
    int GridX;
    int GridY;
    int GridZ;
    float SpaceX;
    float SpaceY;
    float SpaceSlice;
    float SpaceThick;
    short* Voxels;
    // Constructor and Destructor
	  MRIVolData();
	  ~MRIVolData();
};

#endif // MRIVOLDATA_H
