function [S] = adaptiveSettings()
S.totalIters=300;
S.blockSizes = 2;

S.verboseSolver = 0;
S.solver='mosek';
S.finishWithMosek=0;

S.doFullSDPSolutions =4;

S.scsItersVec = 50:50:100;
S.increaseItersVec =50:50:100;
S.increaseItersEvery=5;  % every # 0-level iterations
S.saveWorkspace = false;
S.fileNameComment=[];
S.saveVars =0;
S.warmStart=0;
S.warmStartFile=[];
S.saveVarsInFile = 'testSaveVars'; %no .mat ending needed

S.adaptiveLookback = 3;
S.addItersblockTolerance = @(l) 5e-3.*0.1.^(l>8).*0.1.^(l>16).*0.1.^(l>40) ;%0.01* 10.^(-l/6); %function of current level
S.adaptiveMode = 'down'; 

% scs stopping cirteria
S.eps_abs =1e-8;
S.eps_rel =1e-8;

S.randomInitVars =0;