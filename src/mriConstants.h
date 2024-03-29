#ifndef MRICONST_H
#define MRICONST_H

#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

  // GLOBAL CONSTANTS
  const double kPI = 3.14159;
  // Number Of dimensions
  const int kNumberOfDimensions = 3;
  const int k3DNeighbors = 6;
  // Headers etc.
  const int kHeaderCommentsLines = 13;
  //const int kHeaderCommentsLines = 20;
  const int kInitialMemorySize = 1000;
  // Math
  const double kMathZero = 1.0e-8;
  // File Types
  const int ftPLT = 0;
  const int ftVOL = 1;
  // Adjacency Types
  const int kfacePlusX =  0;
  const int kfaceMinusX = 1;
  const int kfacePlusY =  2;
  const int kfaceMinusY = 3;
  const int kfacePlusZ =  4;
  const int kfaceMinusZ = 5;
  
  // Quantities
  const int kNoQuantity =  -1;
  const int kQtyPositionX =  0;
  const int kQtyPositionY =  1;
  const int kQtyPositionZ =  2;
  const int kQtyConcentration =  3;
  const int kQtyVelocityX =  4;
  const int kQtyVelocityY =  5;
  const int kQtyVelocityZ =  6;
  const int kQtyPressGradientX =  7;
  const int kQtyPressGradientY =  8;
  const int kQtyPressGradientZ =  9;
  const int kQtyPressGradientMod = 10;
  const int kQtyRelPressure =  11;
  const int kQtyVelModule =  12;

  // Threshold Criteria
  const int kCriterionLessThen = 0;
  const int kCriterionGreaterThen = 1;
  const int kCriterionABSLessThen = 2;
  const int kCriterionABSGreaterThen = 3;

  // Directions
  const int kdirX = 0;
  const int kdirY = 1;
  const int kdirZ = 2;
  const int kdirT = 3;
  // Records From Vol Data
  const int kVolAnatomy = 0;
  const int kVolVelocityX = 1;
  const int kVolVelocityY = 2;
  const int kVolVelocityZ = 3;
  const int kPressGradX = 4;
  const int kPressGradY = 5;
  const int kPressGradZ = 6;
  const int kRelPressure = 7;
  
  // Sample Flow Configurations
  const int kZeroVelocity         = 0;
  const int kConstantFlow         = 1;
  const int kPoiseuilleFlow       = 2;
  const int kStagnationFlow       = 3; 
  const int kCylindricalVortex    = 4;
  const int kSphericalVortex      = 5;
  const int kToroidalVortex       = 6;
  const int kTransientFlow        = 7;
  const int kConstantFlowWithStep = 8;
  const int kTaylorVortex         = 9;
  
  // Streamlines Plane
  const int kPlaneXY = 0;
  const int kPlaneYZ = 1;
  const int kPlaneZX = 2;

  // Vortex Criterion
  const int kVortexQ     = 0;
  const int kVortexL2    = 1;
  const int kVortexDelta = 2;

  // Vortex Expansion Thresholding
  const int kHardThresold  = 0;
  const int kSoftThreshold = 1;

  // Smoothing Filter Type
  const int kMedianFilter   = 0;
  const int kMeanFilter     = 1;
  const int kGaussianFilter = 2;

  // Types of inputs
  const int kInputVTK = 0;
  const int kInputPLT = 1;

  // Aternative Typedefs
  typedef const int mriDirection;
  typedef const int mriTemplateType;
  typedef const int mriPlane;
  typedef boost::variate_generator<boost::mt19937, boost::normal_distribution<> > stdRndGenerator;
  
#endif // MRICONST_H
