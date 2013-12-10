#include <math.h>

#include "mriException.h"
#include "mriScan.h"
#include "mriSequence.h"
#include "schMessages.h"

// Make Difference of Scans
// Put the Result Of the Operation in firstScan
void MRISequence::MakeScanDifference(MRIScan &firstScan, const MRIScan &secondScan){
  // set the Tolerance Value
  double DistTol = 1.0e-4;
  // Check Compatibility
  bool scansAreCompatible = (firstScan.cellTotals[0] == secondScan.cellTotals[0])&&
                            (firstScan.cellTotals[1] == secondScan.cellTotals[1])&&
                            (firstScan.cellTotals[2] == secondScan.cellTotals[2]);
  if(!scansAreCompatible) throw MRISequenceException("Scans are not compatible");
  // Write Progress Message
  WriteSchMessage("Computing Scan Difference...");  
  // SubTract Velocity Data
  double diffX,diffY,diffZ;
  for(int loopA=0;loopA<firstScan.totalCellPoints;loopA++){
    // Check They have the same Coordinates
    diffX = fabs((firstScan.cellPoints[loopA].position[0]-secondScan.cellPoints[loopA].position[0])/(firstScan.cellPoints[loopA].position[0]));
    diffY = fabs((firstScan.cellPoints[loopA].position[1]-secondScan.cellPoints[loopA].position[1])/(firstScan.cellPoints[loopA].position[1]));
    diffZ = fabs((firstScan.cellPoints[loopA].position[2]-secondScan.cellPoints[loopA].position[2])/(firstScan.cellPoints[loopA].position[2]));
    // Throw exception if the coordinate is different
    if((diffX>DistTol)||(diffY>DistTol)||(diffZ>DistTol)) throw new MRISequenceException("Error In MakeGlobalDataDifference: Different Coords");
    for(int loopB=0;loopB<3;loopB++){
      firstScan.cellPoints[loopA].velocity[loopB] = firstScan.cellPoints[loopA].velocity[loopB]-secondScan.cellPoints[loopA].velocity[loopB];
    }
  }
}

// Make Average
// Put the Result Of the Operation in firstScan
void MRISequence::MakeScanAverage(int numberOfMeasures, MRIScan &firstScan, const MRIScan &secondScan){
  double DistTol = 1.0e-4;
  // Check Compatibility
  bool scansAreCompatible = (firstScan.cellTotals[0] == secondScan.cellTotals[0])&&
                            (firstScan.cellTotals[1] == secondScan.cellTotals[1])&&
                            (firstScan.cellTotals[2] == secondScan.cellTotals[2]); 
  if(!scansAreCompatible) throw MRISequenceException("Scans are not compatible"); 
  // Write Progress Message
  WriteSchMessage("Computing Scan Average...");
  // SubTract Velocity Data
  double diffX,diffY,diffZ;
  for(int loopA=0;loopA<firstScan.totalCellPoints;loopA++){
    // Check They have the same Coordinates
    diffX = fabs((firstScan.cellPoints[loopA].position[0]-secondScan.cellPoints[loopA].position[0])/(firstScan.cellPoints[loopA].position[0]));
    diffY = fabs((firstScan.cellPoints[loopA].position[1]-secondScan.cellPoints[loopA].position[1])/(firstScan.cellPoints[loopA].position[1]));
    diffZ = fabs((firstScan.cellPoints[loopA].position[2]-secondScan.cellPoints[loopA].position[2])/(firstScan.cellPoints[loopA].position[2]));
    if((diffX>DistTol)||(diffY>DistTol)||(diffZ>DistTol))throw new MRISequenceException("Error In MakeGlobalDataAverage: Different Coords");
    for(int loopB=0;loopB<3;loopB++){
      firstScan.cellPoints[loopA].velocity[loopB] =  
      firstScan.cellPoints[loopA].velocity[loopB]*((double)(numberOfMeasures-1)/(double)numberOfMeasures)+
      secondScan.cellPoints[loopA].velocity[loopB]*(1.0/(double)numberOfMeasures);
    }
  }
}

// Make Difference of Scans
/*void MRISequence::MakeScanDifference(FileName1,FileName2: String;
                           Var GlobalData1: GlobalDataRecord): Integer4;
Var
  ErrorCode: Integer4;
  GlobalData2: GlobalDataRecord;
  Continue: Boolean;
Begin
  {Init Result}
  Result:=MRI_IO_NoError;
  {Open The Two Files}
  {File1}
  ErrorCode:=ReadPltFile(FileName1,TRUE,GlobalData1);
  If (ErrorCode<>MRI_IO_NoError) Then Exit;
  {File2}
  ErrorCode:=ReadPltFile(FileName2,TRUE,GlobalData2);
  If (ErrorCode<>MRI_IO_NoError) Then Exit;
  {Check Compatibility}
  Continue:=(GlobalData1.CellTotals[1]=GlobalData2.CellTotals[1])And
            (GlobalData1.CellTotals[2]=GlobalData2.CellTotals[2])And
            (GlobalData1.CellTotals[3]=GlobalData2.CellTotals[3]);
  If (Continue) Then
  Begin
    {Make Difference: Assign to Global Data 1}
    MakeGlobalDataDifference(GlobalData1,GlobalData2);
  End Else Begin
    MessageDlg('Error: Files are not Compatible.',mtError,[mbOK],0);
    Result:=MRI_DIFF_DataFilesNotCompatible;
    Exit;
  End;
End;

{Read Structure From File}
Function MakePLTAverage(TotalFiles: Integer4;FileNames: StringArray): Integer4;
Var
  LoopA: Integer4;
  ErrorCode: Integer4;
  GlobalData1: GlobalDataRecord;
  GlobalData2: GlobalDataRecord;
  Continue: Boolean;
Begin
  {Init Result}
  Result:=MRI_IO_NoError;

  {Read The First File}
  ErrorCode:=ReadPltFile(FileNames[1],FALSE,GlobalData1);
  If (ErrorCode<>MRI_IO_NoError) Then Exit;

  {Export First Step}
  ErrorCode:=ExportToTECPLOT(ExtractFilePath(FileNames[1])+'\Step_'+IntToStr(1),GlobalData1);
  If (ErrorCode<>MRI_IO_NoError) Then Exit;

  {Loop On All The Other Files}
  For LoopA:=2 To TotalFiles Do
  Begin
    {Read Next File}
    ErrorCode:=ReadPltFile(FileNames[LoopA],FALSE,GlobalData2);
    If (ErrorCode<>MRI_IO_NoError) Then Exit;

    {Check Compatibility}
    Continue:=(GlobalData1.CellTotals[1]=GlobalData2.CellTotals[1])And
              (GlobalData1.CellTotals[2]=GlobalData2.CellTotals[2])And
              (GlobalData1.CellTotals[3]=GlobalData2.CellTotals[3]);
    If (Not(Continue)) Then
    Begin
      MessageDlg('Error: Files are not Compatible.',mtError,[mbOK],0);
      Result:=MRI_AVG_DataFilesNotCompatible;
      Exit;
    End;

    {Make Average}
    MakeGlobalDataAverage(LoopA,GlobalData1,GlobalData2);

    {Print Avergae}
    ErrorCode:=ExportToTECPLOT(ExtractFilePath(FileNames[1])+'\Step_'+IntToStr(LoopA),GlobalData1);
    If (ErrorCode<>MRI_IO_NoError) Then Exit;
  End;
End;

Function WriteGlobalDataDifference(FileName: String;
                                   GlobalDataRef: GlobalDataRecord;
                                   CurrentGlobalData: GlobalDataRecord): Integer4;
Var
  LoopA,LoopB: Integer4;
  MagnitudeDiffFileName: String;
  AngleDoffFileName: String;
  StDevFileName: String;
  MagFile,AngFile,StDevFile: TextFile;
  RefVelocity,CurrentVelocity: Array3Double;
  DiffMagnitude,DiffAngle,StDevValue: Single;
  MagRef,MagCurrent: Single;
  DiffMagnitudeAv,DiffMagnitudeSigma: Single;
  DiffAngleAv,DiffAngleSigma: Single;
  TempValue: Single;
  NumberOfSamples: Integer4;
  CurrentDeviation: Single;
Begin
  {Init Result}
  Result:=MRI_IO_NoError;
  {Find File Names}
  MagnitudeDiffFileName:=FileName+'Mag.dat';
  AngleDoffFileName:=FileName+'Ang.dat';
  StDevFileName:=FileName+'StDev.dat';
  {Assign Files}
  AssignFile(MagFile,MagnitudeDiffFileName);
  AssignFile(AngFile,AngleDoffFileName);
  AssignFile(StDevFile,StDevFileName);
  Rewrite(MagFile);
  Rewrite(AngFile);
  Rewrite(StDevFile);
  {Initialize Statistics}
  NumberOfSamples:=0;
  DiffMagnitudeAv:=0.0;
  DiffMagnitudeSigma:=0.0;
  DiffAngleAv:=0.0;
  DiffAngleSigma:=0.0;
  StDevValue:=0.0;
  {Loop Through All Cells}
  For LoopA:=1 To CurrentGlobalData.TotalCellPoints Do
  Begin
    {Eval Ref Magnitude}
    For LoopB:=1 To 3 Do RefVelocity[LoopB]:=GlobalDataRef.CellPoints[LoopA].Velocity[LoopB];
    {Eval Curren Velocity}
    For LoopB:=1 To 3 Do CurrentVelocity[LoopB]:=CurrentGlobalData.CellPoints[LoopA].Velocity[LoopB];
    CurrentDeviation:=0.0;
    For LoopB:=1 To 3 Do CurrentDeviation:=CurrentDeviation+(1.0/3.0)*(CurrentVelocity[LoopB]-RefVelocity[LoopB]);
    StDevValue:=StDevValue+Sqr(CurrentDeviation);
    {Eval The Average Deviation}
    {Eval Diff Magnitude}
    MagRef:=Sqrt(Sqr(RefVelocity[1])+Sqr(RefVelocity[2])+Sqr(RefVelocity[3]));
    MagCurrent:=Sqrt(Sqr(CurrentVelocity[1])+Sqr(CurrentVelocity[2])+Sqr(CurrentVelocity[3]));
    If ((MagCurrent/GlobalDataRef.MaxVelModule)>5.0e-3) Then
    Begin
      {Increment Number Of Samples}
      Inc(NumberOfSamples);
      {Eval Difference In Magnitude}
      DiffMagnitude:=ABS(MagCurrent-MagRef)/(GlobalDataRef.MaxVelModule);
      {Normalize}
      Normalize3DVector(RefVelocity);
      Normalize3DVector(CurrentVelocity);
      {Eval Angle}
      DiffAngle:=0.0;
      For LoopB:=1 To 3 Do DiffAngle:=DiffAngle+(RefVelocity[LoopB]*CurrentVelocity[LoopB]);
      If (DiffAngle>1.0) Then DiffAngle:=1.0;
      If (DiffAngle<-1.0) Then DiffAngle:=-1.0;
      DiffAngle:=ArcCos(DiffAngle)*(180.0/Pi);
      {Sum Statistics}
      DiffMagnitudeAv:=DiffMagnitudeAv+DiffMagnitude;
      DiffMagnitudeSigma:=DiffMagnitudeSigma+Sqr(DiffMagnitude);
      DiffAngleAv:=DiffAngleAv+DiffAngle;
      DiffAngleSigma:=DiffAngleSigma+Sqr(DiffAngle);
    End;
  End;

  {Eval Average}
  StDevValue:=Sqrt(StDevValue/CurrentGlobalData.TotalCellPoints);

  {Eval Average}
  DiffMagnitudeAv:=(DiffMagnitudeAv/NumberOfSamples);
  DiffMagnitudeSigma:=(DiffMagnitudeSigma/NumberOfSamples);
  DiffAngleAv:=(DiffAngleAv/NumberOfSamples);
  DiffAngleSigma:=(DiffAngleSigma/NumberOfSamples);

  {First and Second Order Statistics}
  TempValue:=DiffMagnitudeSigma-Sqr(DiffMagnitudeAv);
  If (TempValue>MathZero) Then DiffMagnitudeSigma:=Sqrt(TempValue)
  Else DiffMagnitudeSigma:=0.0;

  TempValue:=DiffAngleSigma-Sqr(DiffAngleAv);
  If (TempValue>MathZero) Then DiffAngleSigma:=Sqrt(TempValue)
  Else DiffAngleSigma:=0.0;

  {Write Statistics}
  Writeln(MagFile,Format('%8.4e',[DiffMagnitudeAv]));
  Writeln(MagFile,Format('%8.4e',[DiffMagnitudeSigma]));
  Writeln(AngFile,Format('%8.4e',[DiffAngleAv]));
  Writeln(AngFile,Format('%8.4e',[DiffAngleSigma]));
  Writeln(StDevFile,Format('%8.4e',[StDevValue]));
  {Close File}
  CloseFile(MagFile);
  CloseFile(AngFile);
  CloseFile(StDevFile);
End;

{Read Structure From File}
Function EvalModulusAngleDifference(TotalFiles: Integer4;FileNames: StringArray): Integer4;
Var
  LoopA: Integer4;
  ErrorCode: Integer4;
  GlobalDataRef: GlobalDataRecord;
  CurrentGlobalData: GlobalDataRecord;
  Continue: Boolean;
Begin
  {Init Result}
  Result:=MRI_IO_NoError;

  {Read The First File}
  ErrorCode:=ReadPltFile(FileNames[1],FALSE,GlobalDataRef);
  If (ErrorCode<>MRI_IO_NoError) Then Exit;

  {Loop On All The Other Files}
  For LoopA:=2 To TotalFiles Do
  Begin
    {Read Next File}
    ErrorCode:=ReadPltFile(FileNames[LoopA],FALSE,CurrentGlobalData);
    If (ErrorCode<>MRI_IO_NoError) Then Exit;

    {Check Compatibility}
    Continue:=(GlobalDataRef.CellTotals[1]=CurrentGlobalData.CellTotals[1])And
              (GlobalDataRef.CellTotals[2]=CurrentGlobalData.CellTotals[2])And
              (GlobalDataRef.CellTotals[3]=CurrentGlobalData.CellTotals[3]);
    If (Not(Continue)) Then
    Begin
      MessageDlg('Error: Files are not Compatible.',mtError,[mbOK],0);
      Result:=MRI_AVG_DataFilesNotCompatible;
      Exit;
    End;

    {Make Average}
    WriteGlobalDataDifference('G:\Temp\Data\File_'+IntToStr(LoopA),GlobalDataRef,CurrentGlobalData);
  End;
End;
end.*/
