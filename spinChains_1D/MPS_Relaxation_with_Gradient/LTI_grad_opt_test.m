function [energyFinal,energyInit,Wopt,W0,exitflag,output]=LTI_grad_opt_test(D,N,H,seed) 

 if nargin<4
     seed=1;
 end
 
    maxIters=1000;
 
    dualize=0;
h = H.interactionTerm;
d = H.localDim;

 
k=N-2;

%% Exact diag. to compute isometry
 
HK=sparse(d^(k),d^(k));

for j=1:k-1
    interaction = tensor(speye(d^(j-1)),h,speye(d^(k-j-1)));
    HK = HK + real(interaction);
end

% [eigenvecK,eigenval]=eigs(HK,D,'smallestreal');
[eigenvecK,eigenval]=eig(full(HK));
  
 
yalmip('clear')
 
rho_kp1=sdpvar(d^(k+1));

omega=sdpvar(d^(2)*D);

rng(seed,'twister'); 

doGradOptimization=1;
if D==1
    W0 =2*( rand(1)-0.5)*eigenvecK(:,1)' +2* ( rand(1)-0.5)*eigenvecK(:,2)'      ; 
elseif D==d^k
    W0 = speye(D);
    doGradOptimization=0;
else
    W0 = eigenvecK(:,1:D)' ;
    randomizing = 1;
    if randomizing
        
        warning('RANDOMIZING INITAL MAP')
        [Q,~] = qr(randn(d^k));
        W0 = W0*Q;
    end
end



YalmipOps=sdpsettings('solver','mosek','dualize',dualize , 'verbose',0 );
YalmipOps.scs.max_iters=maxIters;
[energyInit ] =  solveSDPandComputeGrad(k,W0,h,rho_kp1,omega,D,d,YalmipOps);
fprintf('~~~~~ initial relaxation energy from solver : %0.10g \n', energyInit) 

%%

f = @(x) objWrapper(x,k,h,rho_kp1,omega,D,d,YalmipOps);
 
x0 = W0(:); 
%%
 
% options = optimoptions(@fmincon,'Algorithm','interior-point',...
%     'FiniteDifferenceStepSize',ones(numel(x0),1)*1e-6,...
%     'CheckGradients',true,'SpecifyObjectiveGradient',true);
%  [x fval exitflag output] = fmincon(f,x0,[],[],[],[],[],[],[],options);
%%

% dont try to check derivatives if CG dimension is bigger than physical
% dimension. then the derivatives are just noise because the obj should
% never change
if doGradOptimization
    opts  = optimoptions('fminunc','CheckGradients',false,'SpecifyObjectiveGradient',true,'FiniteDifferenceType','central',...
         'FiniteDifferenceStepSize',ones(numel(x0),1)*1e-5,...
         'Display','iter-detailed',...
         'OptimalityTolerance',1e-6,...
         'FunctionTolerance',1e-6 ,...
         'StepTolerance',1e-5,...
         'Algorithm','quasi-newton' );
     % , 'MaxFunctionEvaluations', 20 ... 




    [x,fval,exitflag,output]  = fminunc(f,x0,opts) ;

    energyFinal = -fval;
    Wopt = reshape(x,D,d^k);
else
     x=x0;
     exitflag=[];
     output=[];
     energyFinal = energyInit;
     Wopt = W0;
end

end    
   



function [obj,grad] = objWrapper(x,k,h,rho_kp1,omega,D,d,YalmipOps)
    W = reshape(x,D,d^k); 
    if nargout >1
        [optEnergy,gradEnergy] =   solveSDPandComputeGrad(k,W,h,rho_kp1,omega,D,d,YalmipOps);
        grad= -gradEnergy; % minus sign because we minimize(-f) instead of max(f)
    else
        [optEnergy]      =   solveSDPandComputeGrad(k,W,h,rho_kp1,omega,D,d,YalmipOps);
    end
    obj = -optEnergy;      % minus sign because we minimize(-f) instead of max(f)

end




function [optEnergy,grad,solDiagnost] = solveSDPandComputeGrad(k,W,h,rho_kp1,omega,D,d,YalmipOps)

 
[constraints] = setConstraints(k,rho_kp1,omega,W,D,d);

 
objFunc=trace(kron(h,speye(d^(k-1))) * rho_kp1);
solDiagnost = optimize(constraints,objFunc,YalmipOps);

optEnergy = value( objFunc);

%% compute gradient
if nargout >1

    gamma = struct('L',zeros(d*D),'R',zeros(d*D));
    
    fns = fieldnames(gamma);
    lowerTriangIndices = tril(true(d*D));
    for fni=1:2
        yalmipVersion = yalmip('version');
        switch yalmipVersion
            case {'20210331'}
                G = dual(constraints(fns{fni})) ;
                gamma.(fns{fni}) = (G + G')/2   ;
            case {'20210609'}
                error('for some reason dual var is incorrectly read out in this version =?')
%                 G = zeros(d*D);
%                 G(lowerTriangIndices) = dual(constraints(fns{fni})) ;
%                 gamma.(fns{fni}) = (G + G')/2   ;
        end
    end
    rhoOpt=value(rho_kp1);
    
    IxW=kron(speye(d),W);
    WxI=kron(W,speye(d));
    
    dW = -2* reshape( ...
                ncon({decomposeLegs(gamma.L * WxI * rhoOpt,[D d],d*ones(1,k+1))},{[-1 1, -(k+1):1:-2, 1]} )+ ...
                ncon({decomposeLegs(gamma.R * IxW * rhoOpt,[d D],d*ones(1,k+1))},{[1 -1 1, -(k+1):1:-2]} ), ...
                     D,d^k);
    grad = dW(:);
     
    
    if 0  % test complementarity slackness
       norm((kron(eye(d),gamma.L) + kron(gamma.R,eye(d))) * value(omega))
        
    end
    
    % direct computation of gradient for debugging
%     dW2 = zeros(size(W));
%    
%     for i=1:size(W,1)
%         for j=1:size(W,2)
%             Eij  = circshift([1 zeros(1,size(W,1)-1)]',i-1) * circshift([1 zeros(1,size(W,2)-1)],j-1);
%             IxEij =  kron(speye(d),Eij);
%             EijxI =  kron(Eij,speye(d));
% 
%             dW2(i,j) = -2* trace(...
%                 gamma.L * EijxI * rhoOpt * WxI' + ...
%                 gamma.R * IxEij * rhoOpt * IxW' );
%         end
%     end
% 
%     assert(norm(dW-dW2) < 1e-14)
end
end
 

function [constraints] = setConstraints(k,rho_kp1,omega,W,D,d)
IxW=kron(speye(d),W);
WxI=kron(W,speye(d));

constraints = [...
            trace(rho_kp1) == 1,...
            rho_kp1 >= 0,... 
            omega >= 0,... 
            pTr(rho_kp1,[1],d*ones(1,k+1)) == pTr(rho_kp1,[k+1],d*ones(1,k+1)),...
            [pTr(omega,[1],[d D*d]) == WxI * rho_kp1 * WxI']:'L',...
            [pTr(omega,[2],[d*D d]) == IxW * rho_kp1 * IxW']:'R' ];

end



