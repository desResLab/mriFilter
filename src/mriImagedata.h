#ifndef MRIIMAGEDATA_H
#define MRIIMAGEDATA_H

class mriImageData{
  public:
    // Data Members
    int sizeX;
    int sizeY;
    short* rawData;
    // Constructor and Destructor
    mriImageData();
    ~mriImageData();
};

#endif // MRIIMAGEDATA_H
