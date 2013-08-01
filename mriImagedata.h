#ifndef MRIIMAGEDATA_H
#define MRIIMAGEDATA_H

class MRIImageData
{
public:
  public:
    // Data Members
    int sizeX;
    int sizeY;
    short* rawData;
    // Constructor and Destructor
    MRIImageData();
    ~MRIImageData();
};

#endif // MRIIMAGEDATA_H
