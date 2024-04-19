function [energyFinal,energyInit,MPSopt,mpsInit,output]=MPSgrad_opt_full(D,N,randInitMPS,perturbFinal ) 

if nargin<4
 perturbFinal=0;
end
if nargin <5
    dualize=0;
end
d=2;
 
% heisenberg with sub-lattice rotation
  h = -1* GetTwoSiteH([-1,-1,1,0,0],d);     
% heisenberg
% h =  GetTwoSiteH([-1,-1,-1,0,0],d);     


load('HeisenbergSpinHalf_SubLatticeRotation_Data.mat') 
MPS=MPS{D};

n=N-2;
 
k0=floor(2*log(D)/log(d))+1; % number of spins on which the first coarse-graining map acts. chosen such that d^k0 > D^2

 
if randInitMPS ~= 0
    s=rng(4,'twister');
    for l=1:d
        MPS.AL{l,1}  = MPS.AL{l,1} + (randInitMPS *rand(D));
    
    end
end

[W0 ,L,R,mpsInit]= CGmapsFromMPS_Left(MPS.AL,k0,n);
CGmaps=struct('W0',W0,'mps',mpsInit,'mpsStrings',[]);
CGmaps.L=L;
CGmaps.R=R;
%%

yalmip('clear')
[rho,omega]=declarePrimalVars(n,d,D,k0);

YalmipOps=sdpsettings('solver','mosek','dualize',dualize , 'verbose',0 );
DaP=struct('d',d,'D',D, 'h',h,'N',N,'n',n,'k0',k0); %Dims and Parameters
                  
[energyInit ] =  solveSDPandComputeGrad(DaP,CGmaps,rho,omega,YalmipOps);
fprintf('~~~~~ initial relaxation energy from solver : %0.10g \n', energyInit) 

%%

f = @(x) objWrapper(x,DaP,rho,omega,YalmipOps);
 
x0 = mpsInit(:); 
  
opts  = optimoptions('fminunc','CheckGradients',false,'SpecifyObjectiveGradient',true,'FiniteDifferenceType','central',...
     'FiniteDifferenceStepSize',ones(numel(x0),1)*1e-5,...
     'Display','iter-detailed',...
     'OptimalityTolerance',1e-6,...
     'FunctionTolerance',1e-6 ,...
     'StepTolerance',1e-5,...
     'Algorithm','quasi-newton',...
     'HessUpdate','bfgs'); %bfgs / dfp / steepdesc


[x,fval,exitflag,output]  = fminunc(f,x0,opts) ;

if perturbFinal
   [x,fval,exitflag,output]  = fminunc(f,x + perturbFinal*rand(length(x),1),opts) ;
end


energyFinal = -fval;
MPSopt = reshape(x,[D,d,D]);

end    
   



function [obj,grad] = objWrapper(x, DaP,rho,omega,YalmipOps)
    D=DaP.D;
    d=DaP.d;
    
    n=DaP.n;
    k0=DaP.k0;
    mpsFormat=1;
    
    mps = reshape(x,[D,d,D]); 
    [W0 ,L,R]= CGmapsFromMPS_Left(mps,k0,n,mpsFormat);
     
    CGmaps=struct('W0',W0,'mps',mps,'mpsStrings',[]); 
    CGmaps.L=L;
    CGmaps.R=R;
    
    if nargout >1
        mpsStrings = makeMPSstrings(mps,k0-1); %needed for grad computation
        CGmaps.mpsStrings = mpsStrings;
        [optEnergy,gradEnergy] =   solveSDPandComputeGrad(DaP,CGmaps,rho,omega,YalmipOps);
        grad= -gradEnergy; % minus sign because we minimize(-f) instead of max(f)
    else
        [optEnergy]      =   solveSDPandComputeGrad(DaP,CGmaps,rho,omega,YalmipOps);
    end
    obj = -optEnergy;      % minus sign because we minimize(-f) instead of max(f)

end

function [mpsProd] = contract_mps(mps,k,openInds)
    if nargin<3
        openInds=0;
    end

    % when openInds=1:
    %       -k-1          -k-2
    %        |             |
    %        .-M--M--M--M--.
    %          |  |  |  |
    %         -1 -2 .. -k
    %
    % otherwise reshape into a matrix (bond dims are output)
    
    D=size(mps,1);
    d=size(mps,2);
    assert(size(mps,3)==D,'mps shape should be D,d,D')
    
    indCell={[-(k+1), -1, 1]};
    for l=1:k-2
        indCell = horzcat(indCell, [l, -(l+1), l+1]);
    end
        indCell = horzcat(indCell, [(k-1), -(k) -(k+2)]); % add last entry

    
   if openInds
       mpsProd=ncon(repmat({mps},1,k),indCell);
   else
       mpsProd=ncon(repmat({mps},1,k),indCell,[],[-k-2, -k-1, -k:1:-1]);
       mpsProd=reshape(mpsProd,D^2,d^k); 
   end
end

function [mpsStrings] = makeMPSstrings(mps,m)
        % make mps contractions of all sizes 1..m
        % (keep them as multidimensional arrays (third input to function=1))
        mpsStrings = cell(m,1);
        mpsStrings{1} = permute(mps,[2 1 3]); %permute to have same index convention as for longer strings
        for j=2:m
            mpsStrings{j} = contract_mps(mps,j,1);
        end
end

function [optEnergy,grad,solDiagnost] =  solveSDPandComputeGrad(DaP,CGmaps,rho,omega,YalmipOps)
h=DaP.h;
d=DaP.d;
k0=DaP.k0;
 
[constraints] = setConstraints(DaP,CGmaps,rho,omega);
 
objFunc = trace( kron(h,speye(d^(k0-1))) * rho);
solDiagnost = optimize(constraints,objFunc,YalmipOps);

optEnergy = value( objFunc );

%% compute gradient w.r.t. mps tensor
if nargout >1
    assert(strcmp(yalmip('version'),'20210331'), ...
    'Dual variables convention was checked for YALMIP version 20210331 \n Different convention than in version 20210609. Current version %s',yalmip('version'))

    D=DaP.D;
    d=DaP.d;
    n=DaP.n;
    k0=DaP.k0;
    mpsStrings= CGmaps.mpsStrings;
    
    rhoOpt = value(rho);
    gamma0L = dual(constraints('L0'));
        gamma0L = (gamma0L + gamma0L')/2;
    gamma0R = dual(constraints('R0'));
        gamma0R = (gamma0R + gamma0R')/2;
    
    omegaOpt=cell(n,0);
    gammaL=cell(n,0);
    gammaR=cell(n,0);
    for k=k0:n-1
        omegaOpt{k} = value(omega{k});
        gammaL{k} = dual(constraints(['L',num2str(k)]));
            gammaL{k} = (gammaL{k} + gammaL{k}')/2;
        gammaR{k} = dual(constraints(['R',num2str(k)]));
            gammaR{k} = (gammaR{k} + gammaR{k}')/2;
    end
  
                
    IxW=kron(speye(d),CGmaps.W0);
    WxI=kron(CGmaps.W0,speye(d));
  
    % derivative of first (W0) term (see documentation in
    % 'V0_MPS_grad_opt_test.m'
    Rsandwich = decomposeLegs(gamma0R * IxW * rhoOpt,[d, D, D],d*ones(1,k0+1));
    Rsandwich = ncon({Rsandwich}, {[1, -1, -2, 1, -3:-1:-3-k0+1]});
    
    Lsandwich = decomposeLegs(gamma0L * WxI * rhoOpt,[D, D, d],d*ones(1,k0+1));
    Lsandwich = ncon({Lsandwich}, {[-1, -2, 1, -3:-1:-3-k0+1, 1]});
    
    sandwich = Lsandwich + Rsandwich;  %from here on everything we do to L and R is the same so we can take the sum now.
    
    dM = ncon({sandwich, mpsStrings{k0-1} },{[-1 k0 -2 1:k0-1], [1:k0-1 -3 k0] });

    for pos=2:k0-1
        dM = dM + ncon({sandwich,                          mpsStrings{pos-1},  mpsStrings{k0-pos} },...
                       {[k0, k0+1, 1:pos-1, -2, pos:k0-1], [1:pos-1 k0 -1], [pos:k0-1 -3 k0+1] },...
                       [k0 1:pos-1 pos:k0-1 k0+1]); %contraction order
    end
    
    dM = dM + ncon({sandwich, mpsStrings{k0-1} },{[k0, -3, 1:k0-1, -2], [1:k0-1, k0, -1] });
    
    
    % derivatives of all the next terms
    for k=k0:n-1
        LxI=tensor(CGmaps.L{k},speye(d));
        IxR=tensor(speye(d),CGmaps.R{k});
        
        dM = dM + ncon({decomposeLegs(gammaL{k} * LxI * omegaOpt{k},[D D d],[d D D d]) },...
                       {[-1  1  2, -2 -3  1  2]});
        dM = dM + ncon({decomposeLegs(gammaR{k} * IxR * omegaOpt{k},[d D D],[d D D d]) },...
                       {[ 1  2 -3,  1  2 -1 -2]});
    end
    
    grad= -2* dM(:); % minus sign because of SDP grad formula. factor 2 becuse of d(bra)/d... + d(ket)/d...
    
    
    
     
end
end
 

function [constraints] = setConstraints(DaP,CGmaps,rho,omega)
d=DaP.d;
D=DaP.D;
k0=DaP.k0;
n=DaP.n;
h=DaP.h;
%%
constraints=[];
% positivity
constraints = [constraints, rho  >= 0];
for l=k0:n
    constraints = [constraints, omega{l} >= 0];
end

% local compatability constraint for rho
constraints= [constraints,...
             pTr(rho,1,[d,d^(k0-1),d]) == pTr(rho,3,[d,d^(k0-1),d]) ];
% normalization 
constraints= [constraints, trace(rho) == 1];

% compatibility constraints k=k0
V0xI=tensor(CGmaps.W0,speye(d)); 
IxV0=tensor(speye(d),CGmaps.W0); 
constraints = [constraints,...
                [ pTr(omega{k0},1,[d,D^2,d]) == V0xI * rho * V0xI' ]:['L',num2str(0)],...
                [ pTr(omega{k0},3,[d,D^2,d]) == IxV0 * rho * IxV0' ]:['R',num2str(0)]...
               ];
 
% compatibility constraints k>k0
for k=k0:n-1
    LxI=tensor(CGmaps.L{k},speye(d));
    IxR=tensor(speye(d),CGmaps.R{k});
   constraints = [constraints,...
                 [ pTr(omega{k+1},1,[d,D^2,d]) == LxI * omega{k} * LxI' ]:['L',num2str(k)],...
                 [ pTr(omega{k+1},3,[d,D^2,d]) == IxR * omega{k} * IxR' ]:['R',num2str(k)]...
                 ];
end

end



