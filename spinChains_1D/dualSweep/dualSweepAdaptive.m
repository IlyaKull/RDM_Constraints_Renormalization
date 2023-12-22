function [Es,itersOrder,rigIters,itersWhenLevelsAdded,stats,varsOut] = dualSweepAdaptive(d,n,k0,D,H,sweepOps,sdpOps,CGmaps, CGops,prevSavedVars)

if nargout >5
    outputVars=1;
else 
    outputVars=0;
end

if nargin < 10
    warmStart=0;
    if sweepOps.randomInitVars
        [exprRegister,varsRegister,vars,exprDims] = initVarsRandom(d,n,k0,D,sweepOps.blockSize);
    else
        [exprRegister,varsRegister,vars,exprDims] = initVars(d,n,k0,D,sweepOps.blockSize);
    end
else
    warmStart=1;
    exprRegister = prevSavedVars.exprRegister;
    varsRegister = prevSavedVars.varsRegister;
    vars = prevSavedVars.vars;
    exprDims = prevSavedVars.exprDims;
end

% [vars] = sweepIter(d,n,k0,D,H,CGmaps,exprRegister,varsRegister,vars,exprDims,sdpOps);
% Es=cell2mat(vars.epsil)';
% E1rig = FeasibleDualPt(n,D,k0,H,vars,CGmaps, CGops);

% maximum value 
assert(sweepOps.maxLevel<= (n -k0 +2 -sweepOps.blockSize), 'max allowed value of iterations variable is %d \n',(n -k0 +2 -sweepOps.blockSize) )

solverIn=sdpOps.solver;
 

Es=nan(n-k0+2,sweepOps.totalIters);
rigIters=zeros(1, sweepOps.totalIters);
itersOrder = nan(1,sweepOps.totalIters);
iterblockOrder = nan(1,sweepOps.totalIters);
itersWhenLevelsAdded= zeros(1,sweepOps.totalIters);
stats=struct();
stats.times= Es;
stats.solvers= cell(sweepOps.totalIters,1);
stats.maxIters= zeros(sweepOps.totalIters,1);

if warmStart
    peakLevel=max(prevSavedVars.currentLevels);
    currentIterblock=prevSavedVars.currentIterblock;
    switch sweepOps.adaptiveMode
        case 'updown'
            currentLevels = [0:peakLevel+1, peakLevel:-1:1];
            auxInd=1;
            iShift=-1;
         case 'down'
            currentLevels = [ peakLevel+1:-1:0];
            auxInd=-1;
            iShift=0;
    end
else
    currentIterblock=1;
    switch sweepOps.adaptiveMode
            case 'updown'
                currentLevels = [0 1 ];
                iShift=-1;
            case 'down'
                currentLevels = [ 1 0];
                iShift=0;
    end
    auxInd=1; 
end

for i=1:sweepOps.totalIters
    
    hLevel = currentLevels(mod(auxInd+iShift,numel(currentLevels))+1);
     
       
    itersOrder(i) =hLevel;
    iterblockOrder(i) = currentIterblock;
    
    exprRegister = nan(n,1);
    exprRegister(k0-1:k0-1+hLevel)=0;
    exprRegister(k0-1+hLevel:k0-1+sweepOps.blockSize-1+hLevel)=1;  % one means it is being optimized
    exprRegister(k0-1+sweepOps.blockSize+hLevel:n)=0;              % zero means the expression is frozen
    
    fprintf('~~~~~~~~~~ iteration %d of %d \n',i,sweepOps.totalIters)
    if strncmp(sdpOps.solver,'scs',3)
        fprintf('optimizing expressions: %s, solver: %s, maxIters= %d\n',...
                mat2str(find(exprRegister==1)),sdpOps.solver,sdpOps.scs.max_iters )
    else
        fprintf('optimizing expressions: %s, solver: %s \n',...
                mat2str(find(exprRegister==1)),sdpOps.solver)
    end

    tikker=tic;
    varsRegister = setVarsRegister(exprRegister,varsRegister,n,k0);
    vars = setVars(vars,varsRegister,k0,n,d,D); 
    [vars] = sweepIter(d,n,k0,D,H,CGmaps,exprRegister,varsRegister,vars,exprDims,sdpOps);
    fprintf('result: %0.4g \n', vars.epsil{varsRegister.epsil}')
    Es(:,i)=cell2mat(vars.epsil)';
    stats.times(varsRegister.epsil,i)=toc(tikker);
    stats.solvers{i}=sdpOps.solver;
    stats.maxIters(i)=sdpOps.scs.max_iters;
    
    if hLevel==0 
       
       itersInBlock = itersOrder(iterblockOrder(1:i)==currentIterblock);
       E1sInBlock = Es(1,iterblockOrder(1:i)==currentIterblock);
       E1sLevel0 = E1sInBlock(itersInBlock == 0);
       if numel(E1sLevel0) >= sweepOps.adaptiveLookback
           if sweepOps.finishWithMosek %%%%%%%%%%%%%%%%%%%%%%%%%%
                  sdpOps.solver='mosek'; 
           end                         %%%%%%%%%%%%%%%%%%%%%%%%%%
           E1sLevel0 = E1sLevel0(end-sweepOps.adaptiveLookback+1:end); % look at last sweepOps.adaptiveLookback energies
           peakLevel=max(currentLevels);
           if std(E1sLevel0) < sweepOps.addItersblockTolerance(peakLevel)
               if peakLevel < sweepOps.maxLevel
                   currentIterblock =currentIterblock +1;
                   fprintf('!!!!!!!!!! std of last %d level 0 iters below specified threshold %0.3g. Adding another level: %d \n',...
          sweepOps.adaptiveLookback,sweepOps.addItersblockTolerance(peakLevel),peakLevel+1)
                   switch sweepOps.adaptiveMode
                    case 'updown'
                        currentLevels = [0:peakLevel+1, peakLevel:-1:1];
                        auxInd=1;
                     case 'down'
                        currentLevels = [ peakLevel+1:-1:0];
                        auxInd=-1;
                   end
                   sdpOps.solver=solverIn; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   itersWhenLevelsAdded(i)=peakLevel+1;
               else
                  fprintf('Reached max level specified in sweepOps (%d). Terminating \n',sweepOps.maxLevel) 
                  break
               end
           end
       end
    end
    
    %increase scs max_iters
    if hLevel==0  && ~mod(sum(itersOrder(1:i)==0),sweepOps.increaseItersEvery)
        sdpOps.scs.max_iters=sdpOps.scs.max_iters + sweepOps.increaseIters;
    end
    
     auxInd=auxInd+1;
end


if outputVars
 varsOut=struct();
 varsOut.exprRegister=exprRegister;
 varsOut.varsRegister=varsRegister;
 varsOut.vars=vars;
 varsOut.exprDims=exprDims;
 varsOut.currentLevels=currentLevels;
 varsOut.currentIterblock=currentIterblock;
 varsOut.Es=Es;
 varsOut.itersWhenLevelsAdded=itersWhenLevelsAdded;
end
 
