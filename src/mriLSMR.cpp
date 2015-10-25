/*
Type
  TMatrixTypes = (mtExplicitMat,mtMWMat);
  DivFreeMultTypes = (mtDirect,mtTransposed);
  LSMROptionRecord = Packed Record
    MaxIt: Integer4;
    ANormTol: Double;
    BNormTol: Double;
    WriteMsgs: Boolean;
    IterationPrintStep: Integer4;
  End;

{Solve Least Squares Problem With LSMR}
Function LSMRSolve(LSMROptions: LSMROptionRecord;
                   RowCount: Integer4;
                   ColCount: Integer4;
                   RHS: DoubleArray;
                   Var Converged: Boolean;
                   Var Solution: DoubleArray): String;

implementation

Uses Main_MRI,MRI_MP;

{Application Of the Direct Matrix}
Procedure ApplyDivFreeMatrix(MultType: DivFreeMultTypes;
                             GlobalData: GlobalDataRecord;
                             Vector: DoubleArray;
                             Var Solution: DoubleArray);
Var
  LoopA,LoopB,LoopC,LoopD,LoopE: Integer4;
  ErrorString: String;
  CurrentPos: Integer4;
  TotalSlices,TotalStars: Integer4;
  CurrentIDX,TotalIDX: Integer4;
  TotalStarFaces: Integer4;
  FacesID: Integer4Array;
  FacesCoeffs: DoubleArray;
  CurrentID: Integer4;
Begin
  {Set Current Pos}
  CurrentPos:=0;
  {Loop On The Three Directions}
  For LoopB:=1 To NumberOfDimensions Do
  Begin
    Case LoopB Of
      1: Begin
           {YZ Planes}
           TotalSlices:=GlobalData.CellTotals[1];
           TotalStars:=(GlobalData.CellTotals[2]+1)*(GlobalData.CellTotals[3]+1);
         End;
      2: Begin
           {XZ Planes}
           TotalSlices:=GlobalData.CellTotals[2];
           TotalStars:=(GlobalData.CellTotals[1]+1)*(GlobalData.CellTotals[3]+1);
         End;
      3: Begin
           {XY Planes}
           TotalSlices:=GlobalData.CellTotals[3];
           TotalStars:=(GlobalData.CellTotals[1]+1)*(GlobalData.CellTotals[2]+1);
         End;
    End;
    For LoopC:=1 To TotalSlices Do
    Begin
      CurrentIDX:=(LoopB-1)*TotalSlices+LoopC-1;
      TotalIDX:=NumberOfDimensions*TotalSlices-1;
      If (CurrentIDX Mod (Trunc(TotalIDX/20.0)+1))=0 Then
      Begin
        ShowProgress(Trunc(((CurrentIDX)/(TotalIDX))*100.0));
        Application.ProcessMessages;
      End;
      For LoopD:=1 To TotalStars Do
      Begin
        {Increment Position}
        Inc(CurrentPos);
        {Find Star Shape}
        AssembleStarShape(GlobalData,LoopB{Current Dimension},LoopC{Slice},LoopD{Star},TotalStarFaces,FacesID,FacesCoeffs);
        {Handle Direct and Transposed Multiplication}
        If (MultType=mtDirect) Then
        Begin
          {Direct}
          For LoopE:=1 To TotalStarFaces Do
          Begin
            CurrentID:=FacesID[LoopE];
            Solution[CurrentID]:=Solution[CurrentID]+FacesCoeffs[LoopE]*Vector[CurrentPos];
          End;
        End Else Begin
          {Transposed}
          For LoopE:=1 To TotalStarFaces Do
          Begin
            CurrentID:=FacesID[LoopE];
            try
            Solution[CurrentPos]:=Solution[CurrentPos]+FacesCoeffs[LoopE]*Vector[CurrentID];
            except
              showmessage('Ecco');
            end
          End;
        End;
      End;
    End;
  End;
End;

{Solve Least Squares Problem With LSMR}
Function LSMRSolve(LSMROptions: LSMROptionRecord;
                   RowCount: Integer4;
                   ColCount: Integer4;
                   RHS: DoubleArray;
                   Var Converged: Boolean;
                   Var Solution: DoubleArray): String;
Var
  LoopA: Integer4;
  ErrorString: String;
  {Convergence Check}
  ItCount: Integer4;
  {Variables}
  Alpha,Beta: Double;
  Zeta_Old,ZetaBar,AlphaBar: Double;
  Rho_Old,Rho: Double;
  RhoBar_Old,RhoBar: Double;
  cBar,sBar: Double;
  Betadd,Betad: Double;
  Rhod_Old,TauTilde_Old: Double;
  ThetaTilde,Zeta: Double;
  c,s: Double;
  Theta_New: Double;
  ThetaBar: Double;
  RhoTemp: Double;
  {Variables R}
  BetaHat: Double;
  ThetaTilde_Old,RhoTilde_Old: Double;
  cTilde_Old,sTilde_Old: Double;
  Taud: Double;
  MaxrBar,MinrBar: Double;
  {Vectors}
  uVec,vVec: DoubleArray;
  hVec,hBar,xVec: DoubleArray;
  AvVec,ATuVec: DoubleArray;
  {Norms}
  NormA2,NormA: Double;
  NormAr,CondA: Double;
  NormR,NormX,NormB: Double;
  {Conditions}
  ConvergenceS1,ConvergenceS2,ConvergenceS3: Boolean;
  {Special Cases}
  Num,Den: Double;
Begin
  {Init Result}
  Result:='';

  // ------------------
  // DYNAMIC ALLOCATION
  // ------------------

  SetLength(uVec,RowCount+1);
  SetLength(vVec,ColCount+1);
  SetLength(hVec,ColCount+1);
  SetLength(hBar,ColCOunt+1);
  SetLength(xVec,ColCount+1);
  SetLength(AvVec,RowCount+1);
  SetLength(ATuVec,ColCount+1);

  // -----------------------------------------
  // INITIALIZE VECTORS FROM BIDIAGONALIZATION
  // -----------------------------------------

  {Eval beta_k}
  Beta:=DoEucNorm_1Based(RowCount,RHS);
  {Exit If Beta=0}
  If (ABS(Beta)<MathZero) Then
  Begin
    Result:='Error: Null RHS.';
    Exit;
  End;
  {Eval u_k}
  For LoopA:=1 To RowCount Do uVec[LoopA]:=(RHS[LoopA]/Beta);
  {Apply A^T}

  {Application Of the Direct Matrix}
  ApplyDivFreeMatrix(mtTransposed,GlobalData,uVec,vVec);

  {Eval alpha_k}
  Alpha:=DoEucNorm_1Based(ColCount,vVec);
  {Exit If Alpha=0}
  If (ABS(Beta)<MathZero) Then
  Begin
    Result:='Error: Null Alpha.';
    Exit;
  End;
  {Normalize v_k}
  For LoopA:=1 To ColCount Do vVec[LoopA]:=(vVec[LoopA]/Alpha);

  // ----------------------------------
  // INITIALIZE VARIABLES FOR MAIN LOOP
  // ----------------------------------

  {Eval PsiBar_k}
  ZetaBar:=Alpha*Beta;
  {Eval AlphaBar_k}
  AlphaBar:=Alpha;
  {Eval Rho_km1}
  Rho:=1.0;
  {Eval RhoBar_km1}
  RhoBar:=1.0;
  {Eval CBar_km1}
  cBar:=1.0;
  {Eval SBar_km1}
  sBar:=0.0;
  {Eval h_k}
  For LoopA:=1 To ColCount Do hVec[LoopA]:=vVec[LoopA];
  {Eval hBar_km1}
  For LoopA:=1 To ColCount Do hBar[LoopA]:=0.0;
  {Eval x_k}
  For LoopA:=1 To ColCount Do xVec[LoopA]:=0.0;

  // ------------------------------
  // INITIALIZE VARIABLES FOR ||r||
  // ------------------------------

  Betadd:=Beta;
  Betad:=0.0;
  Rhod_Old:=1.0;
  TauTilde_Old:=0.0;
  ThetaTilde:=0.0;
  Zeta:=0.0;

  // ----------------------
  // Exit if b=0 or A'b = 0
  // ----------------------

  NormAr:=Alpha*Beta;
  If (ABS(NormAr)<MathZero) Then
  Begin
    Result:='Error: RHS has Zero Norm.';
    Exit;
  End;

  {MAIN LOOP}
  NormB:=Beta;
  NormA2:=Sqr(Alpha);
  MaxrBar:=0.0;
  MinrBar:=MaxDouble;
  ItCount:=0;
  Converged:=FALSE;
  While (Not(Converged)And(ItCount<LSMROptions.MaxIt)) Do
  Begin
    {Update Iteration Count}
    Inc(ItCount);

    {Write Progress}
    If LSMROptions.WriteMsgs Then
    Begin
      WriteMsg('Iteration '+IntToStr(ItCount)+'; Residual: '+FloatToStr(NormR));
    End;

    // -----------------
    // BIDIAGONALIZATION
    // -----------------
    {Application Of the Direct Matrix}
    ApplyDivFreeMatrix(mtDirect,GlobalData,vVec,AvVec);


    For LoopA:=1 To RowCount Do uVec[LoopA]:=AvVec[LoopA]-Alpha*uVec[LoopA];
    Beta:=DoEucNorm_1Based(RowCount,uVec);
    If (ABS(Beta)>MathZero) Then
    Begin
      For LoopA:=1 To RowCount Do uVec[LoopA]:=uVec[LoopA]/Beta;

      {Application Of the Transposed Matrix}
      ApplyDivFreeMatrix(mtTransposed,GlobalData,uVec,ATuVec);

      For LoopA:=1 To ColCount Do vVec[LoopA]:=ATuVec[LoopA]-Beta*vVec[LoopA];
      Alpha:=DoEucNorm_1Based(ColCount,vVec);
      If (ABS(Alpha)>MathZero) Then
      Begin
        For LoopA:=1 To ColCount Do vVec[LoopA]:=vVec[LoopA]/Alpha;
      End;
    End;

    // --------------------
    // FIRST PLANE ROTATION
    // --------------------

    Rho_Old:=Rho;
    {Changed AlphaHat With Alpha Bar}
    Rho:=Sqrt(Sqr(AlphaBar)+Sqr(Beta));
    {Changed AlphaHat With Alpha Bar}
    If ABS(Rho)<MathZero Then ShowMessage('Internal: Zero Rho.');
    c:=AlphaBar/Rho;
    s:=Beta/Rho;
    Theta_New:=s*Alpha;
    AlphaBar:=c*Alpha;

    // ---------------------
    // SECOND PLANE ROTATION
    // ---------------------

    RhoBar_Old:=RhoBar;
    Zeta_Old:=Zeta;
    ThetaBar:=sBar*Rho;
    RhoTemp:=cBar*Rho;
    RhoBar:=Sqrt(Sqr(cBar*Rho)+Sqr(Theta_New));
    If ABS(RhoBar)<MathZero Then ShowMessage('Internal: Zero Rho Bar.');
    cBar:=cBar*Rho/RhoBar;
    sBar:=Theta_New/RhoBar;
    Zeta:=cBar*ZetaBar;
    ZetaBar:=-sBar*ZetaBar;

    // ------------------
    // Update h, h_hat, x
    // ------------------

    For LoopA:=1 To ColCount Do hBar[LoopA]:=hVec[LoopA]-(ThetaBar*Rho/(Rho_Old*RhoBar_Old))*hBar[LoopA];
    For LoopA:=1 To ColCount Do xVec[LoopA]:=xVec[LoopA]+(Zeta/(Rho*RhoBar))*hBar[LoopA];
    For LoopA:=1 To ColCount Do hVec[LoopA]:=vVec[LoopA]-(Theta_New/Rho)*hVec[LoopA];

    // --------------
    // Estimate ||r||
    // --------------

    // Apply rotation Q_{k,k+1}.
    BetaHat:=c*Betadd;
    Betadd:=-s*Betadd;

    // Apply rotation Qtilde_{k-1}.
    ThetaTilde_Old:=ThetaTilde;
    RhoTilde_Old:=Sqrt(Sqr(Rhod_Old)+Sqr(ThetaBar));
    If ABS(RhoTilde_Old)<MathZero Then ShowMessage('Internal: Zero RhoTilde_Old.');
    cTilde_Old:=Rhod_Old/RhoTilde_Old;
    sTilde_Old:=ThetaBar/RhoTilde_Old;
    ThetaTilde:=sTilde_Old*RhoBar;
    Rhod_Old:=cTilde_Old*RhoBar;
    Betad:=-sTilde_Old*Betad+cTilde_Old*BetaHat;

    // Forward Substitution
    TauTilde_Old:=(Zeta_Old-ThetaTilde_Old*TauTilde_Old)/RhoTilde_Old;
    If ABS(Rhod_Old)<MathZero Then ShowMessage('Internal: Zero Rhod_Old.');
    Taud:=(Zeta-ThetaTilde*TauTilde_Old)/Rhod_Old;
    NormR:=Sqrt(Sqr(Betad-Taud)+Sqr(Betadd));

    // Estimate ||A||.
    NormA2:=NormA2+Sqr(Beta);
    NormA:=Sqrt(normA2);
    NormA2:=NormA2+Sqr(Alpha);

    // Estimate cond(A).
    MaxrBar:=max(MaxrBar,RhoBar_Old);
    If (ItCount>1) Then MinrBar:=min(MinrBar,RhoBar_Old);
    CondA:=max(MaxrBar,RhoTemp)/min(MinrBar,RhoTemp);

    // Compute norms for convergence testing.
    NormAr:=ABS(ZetaBar);
    NormX:=DoEucNorm_1Based(ColCount,xVec);

    {Check Convergence}
    ConvergenceS1:=(NormR<=(LSMROptions.BNormTol*NormB+LSMROptions.ANormTol*NormA*NormX));
    ConvergenceS2:=(NormAr<=(LSMROptions.ANormTol*NormA*NormR));
    ConvergenceS3:=(NormAr<=(LSMROptions.ANormTol*NormA*NormR));
    Converged:=(ConvergenceS1)Or(ConvergenceS2)Or(ConvergenceS3);

  End;

  {Copy Solution Vector If Converged}
  SetLength(Solution,ColCount+1);
  For LoopA:=1 To ColCount Do Solution[LoopA]:=xVec[LoopA];

  // ------------
  // DEALLOCATION
  // ------------

  SetLength(uVec,0);
  SetLength(vVec,0);
  SetLength(hVec,0);
  SetLength(hBar,0);
  SetLength(xVec,0);
  SetLength(AvVec,0);
  SetLength(ATuVec,0);
  FreeMemory(uVec);
  FreeMemory(vVec);
  FreeMemory(hVec);
  FreeMemory(hBar);
  FreeMemory(xVec);
  FreeMemory(AvVec);
  FreeMemory(ATuVec);
End;

end.*/
