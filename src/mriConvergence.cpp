#include <math.h>
#include "mriUtils.h"
#include "mriCoordItem.cpp"



// FIND MAG AND ANGLE ERROR DISTRIBUTION FOR PERBURBING COORDINATES
void FindDistributions(int vecSize, double* RefVector, MRICoordItem* PerturbedVector, double* &MagDistribution, double* &AngleDistribution){
  // Allocate
  MagDistribution = new double[vecSize];
  AngleDistribution = new double[vecSize];
	double currentVec[3] = {0.0};
  // Find Distribution of Magnitudes and Angles
	double currentMag,refMag,currentAngle;
	double intProd;
  for(int loopA=0;loopA<vecSize;loopA++){
    currentMag = sqrt((PerturbedVector[loopA].x * PerturbedVector[loopA].x)+
		                  (PerturbedVector[loopA].y * PerturbedVector[loopA].y)+
											(PerturbedVector[loopA].z * PerturbedVector[loopA].z));
    refMag = MRIUtils::Do3DEucNorm(RefVector);
		// Store Current Vector
    currentVec[0] = PerturbedVector[loopA].x;
    currentVec[1] = PerturbedVector[loopA].y;
    currentVec[2] = PerturbedVector[loopA].z;
		// Normalize Vector
    MRIUtils::Normalize3DVector(RefVector);
    MRIUtils::Normalize3DVector(currentVec);
    intProd = currentVec[0]*RefVector[0]+currentVec[1]*RefVector[1]+currentVec[2]*RefVector[2];
    if(intProd>1.0) intProd = 1.0;
    if(intProd<-1.0) intProd = -1.0;
    currentAngle = acos(intProd) * (180.0/kPI);
    MagDistribution[loopA] = currentMag - refMag;
    AngleDistribution[loopA] = currentAngle;		
	}
  // Eval Average Values and Standard Deviations
  double MagAverage = 0.0;
  double MagStDev = 0.0;
  double AngleAv = 0.0;
  double AngleStDev = 0.0;
  for(int loopA=0;loopA<vecSize;loopA++){
    MagAverage = MagAverage + MagDistribution[loopA];
    MagStDev =   MagStDev + (MagDistribution[loopA]*MagDistribution[loopA]);
    AngleAv =    AngleAv + AngleDistribution[loopA];
    AngleStDev = AngleStDev + (AngleDistribution[loopA]*AngleDistribution[loopA]);		
	}
  // Eval St Dev
  MagStDev = sqrt(MagStDev - (MagAverage*MagAverage));
  AngleStDev = sqrt(AngleStDev - (AngleAv*AngleAv));
  // Renormalize The Distributions
  for(int loopA=0;loopA<vecSize;loopA++){
    //MagDistribution[LoopA]:=((MagDistribution[LoopA]-MagAverage)/MagStDev);
    //AngleDistribution[LoopA]:=((AngleDistribution[LoopA]-AngleAv)/AngleStDev);
    MagDistribution[loopA] = ((MagDistribution[loopA]-MagAverage));
    AngleDistribution[loopA] = ((AngleDistribution[loopA]-AngleAv));		
	} 
}

/*
// Eval Polar Average
Function EvalPolarAverage(FirstVector,SecondVector: Array3Double): Array3Double;
Var
  FirstMod,SecondMod,AvMod: Double;
  AvVector: Array3Double;
Begin
  {Find the Modulus}
  FirstMod:=DoEucNorm(FirstVector);
  SecondMod:=DoEucNorm(SecondVector);
  AvMod:=0.5*(FirstMod+SecondMod);

  {Normalize The Two Vectors}
  Normalize3DVector(FirstVector);
  Normalize3DVector(SecondVector);

  {Eval Average Unit Vector}
  Result[1]:=0.5*(FirstVector[1]+SecondVector[1]);
  Result[2]:=0.5*(FirstVector[2]+SecondVector[2]);
  Result[3]:=0.5*(FirstVector[3]+SecondVector[3]);

  {Normalize}
  Normalize3DVector(Result);

  {Multiply by the Modulus}
  Result[1]:=Result[1]*AvMod;
  Result[2]:=Result[2]*AvMod;
  Result[3]:=Result[3]*AvMod;
End;

{Find Average Of two Vectors}
Function FindAverageOfTwoVectors(ConvergenceType: StatConvTypes;
                                 FirstVector,SecondVector: Array3Double): Array3Double;
Begin
  If (ConvergenceType=ccComponents) Then
  Begin
    {Components}
    Result[1]:=0.5*(FirstVector[1]+SecondVector[1]);
    Result[2]:=0.5*(FirstVector[2]+SecondVector[2]);
    Result[3]:=0.5*(FirstVector[3]+SecondVector[3]);
  End Else Begin
    {Polar}
    Result:=EvalPolarAverage(FirstVector,SecondVector);
  End;
End;

{Find Average Vector}
Function FindAverageVector(ConvergenceType: StatConvTypes;
                           NumberOfPerturbation: Integer4;
                           PerturbedVectors: CoordItemArray): Array3Double;
Var
  LoopA,LoopB: Integer4;
  CurrentVec: Array3Double;
  AvMag: Double;
Begin
  If (ConvergenceType=ccComponents) Then
  Begin
    {Components}
    Result[1]:=0.0;
    Result[2]:=0.0;
    Result[3]:=0.0;
    {Start the Averaging Process}
    For LoopA:=1 To NumberOfPerturbation Do
    Begin
      Result[1]:=Result[1]+PerturbedVectors[LoopA].X/NumberOfPerturbation;
      Result[2]:=Result[2]+PerturbedVectors[LoopA].Y/NumberOfPerturbation;
      Result[3]:=Result[3]+PerturbedVectors[LoopA].Z/NumberOfPerturbation;
    End;
  End Else Begin
    {Polar}
    Result[1]:=0.0;
    Result[2]:=0.0;
    Result[3]:=0.0;
    AvMag:=0.0;
    {Start the Averaging Process}
    For LoopA:=1 To NumberOfPerturbation Do
    Begin
      CurrentVec[1]:=PerturbedVectors[LoopA].X;
      CurrentVec[2]:=PerturbedVectors[LoopA].Y;
      CurrentVec[3]:=PerturbedVectors[LoopA].Z;
      AvMag:=AvMag+DoEucNorm(CurrentVec)/NumberOfPerturbation;
      Normalize3DVector(CurrentVec);
      For LoopB:=1 To 3 Do Result[LoopB]:=Result[LoopB]+CurrentVec[LoopB];
    End;
    {Write Result}
    Normalize3DVector(Result);
    For LoopA:=1 To 3 Do Result[LoopA]:=Result[LoopA]*AvMag;
  End;
End;

{Eval Angle}
Function EvalAngleDifference(GlobalAvVector,CurrentAvVector: Array3Double): Double;
Var
  InnerProd: Double;
Begin
  {Normalize}
  Normalize3DVector(GlobalAvVector);
  Normalize3DVector(CurrentAvVector);
  InnerProd:=(GlobalAvVector[1]*CurrentAvVector[1]+
              GlobalAvVector[2]*CurrentAvVector[2]+
              GlobalAvVector[3]*CurrentAvVector[3]);
  If InnerProd>1.0 Then InnerProd:=1.0;
  If InnerProd<-1.0 Then InnerProd:=-1.0;
  Result:=ArcCos(InnerProd)*(180.0/Pi);
End;

{Eval Difference In Module}
Function EvalModuleDifference(RefVector,OtherVector: Array3Double): Double;
Begin
  Result:=DoEucNorm(OtherVector)-DoEucNorm(RefVector);
End;

{Evaluate Convergence}
Procedure EvalConvergence(ConvergenceType: StatConvTypes;
                          StartingVector: Array3Double;
                          NumberOfPerturbation: Integer4;
                          PerturbedVectors: CoordItemArray;
                          Var ConvMag: DoubleArray;
                          Var ConvAngle: DoubleArray);
Var
  LoopA,LoopB: Integer4;
  GlobalAvVector: Array3Double;
  FirstVector,SecondVector: Array3Double;
  CurrentAvVector: Array3Double;
Begin
  {Allocate}
  SetLength(ConvMag,NumberOfPerturbation+1);
  SetLength(ConvAngle,NumberOfPerturbation+1);
  {Eval Partial Average}
  For LoopA:=1 To NumberOfPerturbation Do
  Begin
    {Init}
    CurrentAvVector:=FindAverageVector(ConvergenceType,LoopA,PerturbedVectors);
    {Eval Difference in Module}
    ConvMag[LoopA]:=ABS(EvalModuleDifference(StartingVector,CurrentAvVector));
    {Eval Difference in Angle}
    ConvAngle[LoopA]:=ABS(EvalAngleDifference(StartingVector,CurrentAvVector));
  End;
End;

Procedure PrintArrayToFile(FileName: String;
                           Size: Integer4;
                           Vector: DoubleArray);
Var
  LoopA: Integer4;
  OutFile: TextFile;
Begin
  AssignFile(OutFile,FileName);
  Rewrite(OutFile);
  For LoopA:=1 To Size Do
  Begin
    Writeln(OutFile,Format('%8.4e',[Vector[LoopA]]));
  End;
  CloseFile(OutFile);
End;

{Main Routine}
Function PerformConvergenceTest: Integer4;
Var
  LoopA,LoopB: Integer4;
  NumberOfPerturbation: Integer4;
  MonteCarloRuns: Integer4;
  StartingVector: Array3Double;
  State: HqrndState;
  PerturbedVectors: CoordItemArray;
  MagDistribution,AngleDistribution: DoubleArray;
  {Convergence}
  ConvMag_1,ConvAngle_1: DoubleArray;
  ConvMag_2,ConvAngle_2: DoubleArray;
  {Average Values}
  AvMagDistribution: DoubleArray;
  AvAngleDistribution: DoubleArray;
  AvConvMag_1: DoubleArray;
  AvConvAngle_1: DoubleArray;
  AvConvMag_2: DoubleArray;
  AvConvAngle_2: DoubleArray;
Begin
  {Init Result}
  Result:=MRI_IO_NoError;

  {Set Parameters}
  NumberOfPerturbation:=500;
  MonteCarloRuns:=1000;

  {Initialize Random Generator}
  //HQRNDSeed(1,1,State);
  HQRNDRandomize(State);

  {Allocate and Initialize}
  SetLength(AvMagDistribution,NumberOfPerturbation+1);
  For LoopA:=1 To NumberOfPerturbation Do AvMagDistribution[LoopA]:=0.0;
  SetLength(AvAngleDistribution,NumberOfPerturbation+1);
  For LoopA:=1 To NumberOfPerturbation Do AvAngleDistribution[LoopA]:=0.0;
  SetLength(AvConvMag_1,NumberOfPerturbation+1);
  For LoopA:=1 To NumberOfPerturbation Do AvConvMag_1[LoopA]:=0.0;
  SetLength(AvConvAngle_1,NumberOfPerturbation+1);
  For LoopA:=1 To NumberOfPerturbation Do AvConvAngle_1[LoopA]:=0.0;
  SetLength(AvConvMag_2,NumberOfPerturbation+1);
  For LoopA:=1 To NumberOfPerturbation Do AvConvMag_2[LoopA]:=0.0;
  SetLength(AvConvAngle_2,NumberOfPerturbation+1);
  For LoopA:=1 To NumberOfPerturbation Do AvConvAngle_2[LoopA]:=0.0;

  {Run MonteCarlo Tests}
  For LoopA:=1 To MonteCarloRuns Do
  Begin
    {Select a Vector}
    StartingVector:=GenerateRandomVector(State);

    {Perturb its coordinates with Gaussian Noise}
    PerturbVectorGaussian(State,StartingVector,NumberOfPerturbation,PerturbedVectors);

    {Find the Distribution of Magnitudes and angles}
    FindDistributions(StartingVector,NumberOfPerturbation,PerturbedVectors,MagDistribution,AngleDistribution);

    {Make The Average and Eval Convergence}
    EvalConvergence(ccComponents,StartingVector,NumberOfPerturbation,PerturbedVectors,ConvMag_1,ConvAngle_1);
    EvalConvergence(ccPolar,StartingVector,NumberOfPerturbation,PerturbedVectors,ConvMag_2,ConvAngle_2);

    {Sum to Average Values}
    For LoopB:=1 To NumberOfPerturbation Do
    Begin
      AvMagDistribution[LoopB]:=AvMagDistribution[LoopB]+MagDistribution[LoopB];
      AvAngleDistribution[LoopB]:=AvAngleDistribution[LoopB]+AngleDistribution[LoopB];
      AvConvMag_1[LoopB]:=AvConvMag_1[LoopB]+ConvMag_1[LoopB];
      AvConvAngle_1[LoopB]:=AvConvAngle_1[LoopB]+ConvAngle_1[LoopB];
      AvConvMag_2[LoopB]:=AvConvMag_2[LoopB]+ConvMag_2[LoopB];
      AvConvAngle_2[LoopB]:=AvConvAngle_2[LoopB]+ConvAngle_2[LoopB];
    End;
  End;
  {Divide by the Number of Runs}
  For LoopB:=1 To NumberOfPerturbation Do
  Begin
    AvMagDistribution[LoopB]:=(AvMagDistribution[LoopB]/MonteCarloRuns);
    AvAngleDistribution[LoopB]:=(AvAngleDistribution[LoopB]/MonteCarloRuns);
    AvConvMag_1[LoopB]:=(AvConvMag_1[LoopB]/MonteCarloRuns);
    AvConvAngle_1[LoopB]:=(AvConvAngle_1[LoopB]/MonteCarloRuns);
    AvConvMag_2[LoopB]:=(AvConvMag_2[LoopB]/MonteCarloRuns);
    AvConvAngle_2[LoopB]:=(AvConvAngle_2[LoopB]/MonteCarloRuns);
  End;
  {Print to Files}
  PrintArrayToFile('G:\Temp\Data\AvMagDistribution.dat',NumberOfPerturbation,AvMagDistribution);
  PrintArrayToFile('G:\Temp\Data\AvAngleDistribution.dat',NumberOfPerturbation,AvAngleDistribution);
  PrintArrayToFile('G:\Temp\Data\AvConvMag_1.dat',NumberOfPerturbation,AvConvMag_1);
  PrintArrayToFile('G:\Temp\Data\AvConvAngle_1.dat',NumberOfPerturbation,AvConvAngle_1);
  PrintArrayToFile('G:\Temp\Data\AvConvMag_2.dat',NumberOfPerturbation,AvConvMag_2);
  PrintArrayToFile('G:\Temp\Data\AvConvAngle_2.dat',NumberOfPerturbation,AvConvAngle_2);
End;

{Perform Test 01}
Function PerformTest01: Integer4;
Var
  NumberOfPerturbation: Integer4;
  StartingVector: Array3Double;
  State: hqrndState;
  PerturbedVectors: CoordItemArray;
  ConvMag,ConvAngle: DoubleArray;
Begin
  {Init Result}
  Result:=MRI_IO_NoError;

  {Set Params}
  NumberOfPerturbation:=1000;

  {Set Vector}
  StartingVector[1]:=0.0;
  StartingVector[2]:=0.0;
  StartingVector[3]:=1.0;

  {Randomize}
  HQRNDRandomize(State);

  {Perturb its coordinates with Gaussian Noise}
  PerturbVectorGaussian(State,StartingVector,NumberOfPerturbation,PerturbedVectors);

  {Make The Average and Eval Convergence}
  EvalConvergence(ccComponents,StartingVector,NumberOfPerturbation,PerturbedVectors,ConvMag,ConvAngle);

  {Print Data}
  PrintArrayToFile('G:\Temp\Data\ConvMag.dat',NumberOfPerturbation,ConvMag);
  PrintArrayToFile('G:\Temp\Data\ConvAngle.dat',NumberOfPerturbation,ConvAngle);
End;

{Perform Test 02}
Function PerformTest02: Integer4;
Var
  LoopA: Integer4;
  MonteCarloRuns: Integer4;
  Average: Double;
  State: HqrndState;
Begin
  {Randomize}
  HQRNDRandomize(State);

  {Init The MC Runs}
  MonteCarloRuns:=100000;

  {Set the Average}
  Average:=0.0;
  For LoopA:=1 To MonteCarloRuns Do
  Begin
    Average:=Average+(Sqrt(Sqr(HQRNDNormal(State)*0.1)+
                           Sqr(HQRNDNormal(State)*0.1)+
                           Sqr(HQRNDNormal(State)*0.1)))/MonteCarloRuns;
  End;
  ShowMessage('Average: '+FloatToStr(Average));

End;

end.*/
