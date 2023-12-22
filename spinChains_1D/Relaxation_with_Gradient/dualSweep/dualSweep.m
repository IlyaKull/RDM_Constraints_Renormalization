function [Es,itersOrder,rigIters,vars] = dualSweep(d,n,k0,D,H,sweepOps,sdpOps,CGmaps, CGops)

[exprRegister,varsRegister,vars,exprDims] = initVars(d,n,k0,D,sweepOps.blockSize);

% [vars] = sweepIter(d,n,k0,D,H,CGmaps,exprRegister,varsRegister,vars,exprDims,sdpOps);
% Es=cell2mat(vars.epsil)';
% E1rig = FeasibleDualPt(n,D,k0,H,vars,CGmaps, CGops);

maxLevel =sweepOps.maxLevel;
oneSweepCycle=[];
switch sweepOps.cycleType
    case 'linear'
        oneSweepCycle=[oneSweepCycle,1:maxLevel, maxLevel-1:-1:0]; 
    case 'yoyo1'
        for ccl= 1:maxLevel  
            oneSweepCycle = [oneSweepCycle,1:ccl,ccl-1:-1:0];
        end
    case 'yoyo2'
        for ccl= 1:maxLevel  
            oneSweepCycle = [oneSweepCycle,1:ccl,ccl-1:-1:0];
            oneSweepCycle = [oneSweepCycle,1:ccl,ccl-1:-1:0];
        end
    case 'yoyo3'
        for ccl= 1:maxLevel  
            oneSweepCycle = [oneSweepCycle,1:ccl,ccl-1:-1:0];
            oneSweepCycle = [oneSweepCycle,1:ccl,ccl-1:-1:0];
            oneSweepCycle = [oneSweepCycle,1:ccl,ccl-1:-1:0];
        end
case 'yoyo4'
        for ccl= 1:maxLevel  
            oneSweepCycle = [oneSweepCycle,1:ccl,ccl-1:-1:0];
            oneSweepCycle = [oneSweepCycle,1:ccl,ccl-1:-1:0];
            oneSweepCycle = [oneSweepCycle,1:ccl,ccl-1:-1:0];
            oneSweepCycle = [oneSweepCycle,1:ccl,ccl-1:-1:0];
        end
case 'yoyo5'
        for ccl= 1:maxLevel  
            oneSweepCycle = [oneSweepCycle,1:ccl,ccl-1:-1:0];
            oneSweepCycle = [oneSweepCycle,1:ccl,ccl-1:-1:0];
            oneSweepCycle = [oneSweepCycle,1:ccl,ccl-1:-1:0];
            oneSweepCycle = [oneSweepCycle,1:ccl,ccl-1:-1:0];
            oneSweepCycle = [oneSweepCycle,1:ccl,ccl-1:-1:0];
        end
    case 'oneby1'
        for ccl= 1:maxLevel  
            oneSweepCycle = [oneSweepCycle,ccl:-1:0];
        end
    case 'oneby2'
        for ccl= 1:maxLevel  
            oneSweepCycle = [oneSweepCycle,ccl:-1:0];
            oneSweepCycle = [oneSweepCycle,ccl:-1:0];
        end
    case 'oneby3'
        for ccl= 1:maxLevel  
            oneSweepCycle = [oneSweepCycle,ccl:-1:0];
            oneSweepCycle = [oneSweepCycle,ccl:-1:0];
            oneSweepCycle = [oneSweepCycle,ccl:-1:0];       
        end
 case 'oneby4'
        for ccl= 1:maxLevel  
            oneSweepCycle = [oneSweepCycle,ccl:-1:0];
            oneSweepCycle = [oneSweepCycle,ccl:-1:0];
            oneSweepCycle = [oneSweepCycle,ccl:-1:0];       
            oneSweepCycle = [oneSweepCycle,ccl:-1:0];       
        end
 case 'oneby5'
        for ccl= 1:maxLevel  
            oneSweepCycle = [oneSweepCycle,ccl:-1:0];
            oneSweepCycle = [oneSweepCycle,ccl:-1:0];
            oneSweepCycle = [oneSweepCycle,ccl:-1:0];       
            oneSweepCycle = [oneSweepCycle,ccl:-1:0];
            oneSweepCycle = [oneSweepCycle,ccl:-1:0];       
        end
    otherwise
        error('cycle should be either *linear*, *oneby1...5*, or *yoyo1...5* ')
end


if numel(oneSweepCycle) >= sweepOps.totalIters
    sweepOps.totalIters = numel(oneSweepCycle)+1;
end

itersOrder=[0];
while numel(itersOrder) < sweepOps.totalIters
    itersOrder=[itersOrder,oneSweepCycle];
end
% the last iteration should do the 0th levell
finalIter=find(itersOrder(sweepOps.totalIters:end)==0,1) + sweepOps.totalIters - 1;
itersOrder=itersOrder(1:finalIter);
assert(itersOrder(1)==0,'iterations do not start from 0th level')
assert(itersOrder(end)==0,'iterations do not end with  0th level')
% itersOrder=repmat(oneSweepCycle,1,sweepOps.num_sweeps);

% maximum value for II is 
assert(max(itersOrder)<= (n -k0 +2 -sweepOps.blockSize), 'max allowed value of iterations variable is %d \n',(n -k0 +2 -sweepOps.blockSize) )

Es=nan(n-k0+2,numel(itersOrder));
rigIters=zeros(1, numel(itersOrder));

stopCritMet=0;
lastIter=0;

for i=1:numel(itersOrder)
     hLevl= itersOrder(i);
%     if lastIter || i==numel(iterationsOrder) %solve with mosek on last iter
%         hLevl = 0;
%         sdpOps.solver='mosek';
%     else
%         hLevl= iterationsOrder(i); 
%     end
%     
    if i==sweepOps.switchToMosekAfter
        sdpOps.solver='mosek';
    end
    
    
    
    exprRegister = nan(n,1);
    exprRegister(k0-1:k0-1+hLevl)=0;
    exprRegister(k0-1+hLevl:k0-1+sweepOps.blockSize-1+hLevl)=1;  % one means it is being optimized
    exprRegister(k0-1+sweepOps.blockSize+hLevl:n)=0;              % zero means the expression is frozen
    
    fprintf('~~~~~~~~~~ iteration %d of %d \n',i,numel(itersOrder))
    if strncmp(sdpOps.solver,'scs',3)
        fprintf('optimizing expressions: %s, solver: %s, maxIters= %d\n',...
                mat2str(find(exprRegister==1)),sdpOps.solver,sdpOps.scs.max_iters )
    else
        fprintf('optimizing expressions: %s, solver: %s \n',...
                mat2str(find(exprRegister==1)),sdpOps.solver)
    end

    varsRegister = setVarsRegister(exprRegister,varsRegister,n,k0);
    vars = setVars(vars,varsRegister,k0,n,d,D); 
    [vars] = sweepIter(d,n,k0,D,H,CGmaps,exprRegister,varsRegister,vars,exprDims,sdpOps);
    fprintf('result: %0.4g \n', vars.epsil{varsRegister.epsil}')
    Es(:,i)=cell2mat(vars.epsil)';
    
    if (sum(itersOrder(1:i)==0) > sweepOps.stopTolRange ) && ~ lastIter
        e1s = Es(1,1:i);
        e1s = e1s(itersOrder(1:i)==0);
        e1sSTDdev = std(e1s(end-sweepOps.stopTolRange+1:end)); 
        stopCritMet = e1sSTDdev <= sweepOps.stopTol  ;
        if stopCritMet
            fprintf('Standard deviation of last %d solutions of level 0 is %0.2g which is less than specified tol (%0.2g) \n',...
                sweepOps.stopTolRange,e1sSTDdev,sweepOps.stopTol )
            fprintf('************ Terminating after next iteration \n')
        end
    end

    
      makeFeasible =... % conditions for forcing the variables to be strictly feasible
                    ...     % at level 0 every sweepOps.forceFeasabilityEvery times
        (itersOrder(i)==0 && ~mod(sum(itersOrder(1:i)==0),sweepOps.forceFeasabilityEvery) ) ||... 
                    ...     % or at last iteration
        lastIter || i==numel(itersOrder);
    if makeFeasible
        [E1rig,vars] =  FeasibleDualPt(n,D,k0,H,vars,CGmaps, CGops) ;
        fprintf('~~~~~~~~~~~~ verifying rigorous solution:\n difference is: %0.4g \n', E1rig-Es(1,i))
        Es(1,i)=E1rig;
        rigIters(i)=1;
    end     
    
    %increase scs max_iters
    if itersOrder(i)==0 && ~mod(sum(itersOrder(1:i)==0),sweepOps.increaseItersEvery)
        sdpOps.scs.max_iters=sdpOps.scs.max_iters + sweepOps.increaseIters;
    end

     if lastIter
         itersOrder = itersOrder(1:i);
         break
     end
     if stopCritMet
         lastIter=1 ;
     end
end
 
