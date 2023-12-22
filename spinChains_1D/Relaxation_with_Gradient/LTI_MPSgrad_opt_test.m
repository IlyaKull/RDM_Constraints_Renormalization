function [energyFinal,energyInit,MPSopt,mpsInit]=LTI_MPSgrad_opt_test(bondDim,N,randomizeMPS,dualize) 

 
if nargin <4
    dualize=0;
end
d=2;
 
% heisenberg with sub-lattice rotation
% h = -1* GetTwoSiteH([-1,-1,1,0,0],d);     
% heisenberg
h =  GetTwoSiteH([-1,-1,-1,0,0],d);     


load('HeisenbergSpinHalf_SubLatticeRotation_Data.mat') 
MPS=MPS{bondDim};

k=N-2;
 
yalmip('clear')
 
rho_kp1=sdpvar(d^(k+1));

D = bondDim^2;
omega=sdpvar(d^(2)*D);

 

mpsInit=zeros(bondDim,d,bondDim);



if randomizeMPS ~= 0
    s=rng(4,'twister')
    for l=1:d
        mpsInit(:,l,:) = MPS.AL{l,1} + (randomizeMPS *rand(bondDim));
    
    end
else
    for l=1:d
        mpsInit(:,l,:)=MPS.AL{l,1};
    end
    
end

 

W0 = contract_mps(mpsInit,k);
 
YalmipOps=sdpsettings('solver','mosek','dualize',dualize , 'verbose',0 );
 
[energyInit ] =  solveSDPandComputeGrad(k,W0,h,rho_kp1,omega,D,d,mpsInit,YalmipOps);
fprintf('~~~~~ initial relaxation energy from solver : %0.10g \n', energyInit) 

%%

f = @(x) objWrapper(x,k,h,rho_kp1,omega,D,d,YalmipOps);
 
x0 = mpsInit(:); 
%%
 
% options = optimoptions(@fmincon,'Algorithm','interior-point',...
%     'FiniteDifferenceStepSize',ones(numel(x0),1)*1e-6,...
%     'CheckGradients',true,'SpecifyObjectiveGradient',true);
%  [x fval exitflag output] = fmincon(f,x0,[],[],[],[],[],[],[],options);
%%
% dont try to check derivatives if CG dimension is bigger than physical
% dimension. then the derivatives are just noise because the obj should
% never change

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
MPSopt = reshape(x,[bondDim,d,bondDim]);

end    
   



function [obj,grad] = objWrapper(x,k,h,rho_kp1,omega,D,d,YalmipOps)
    bondDim=sqrt(D);
    mps = reshape(x,[bondDim,d,bondDim]); 
    W = contract_mps(mps,k);
    
    if nargout >1
        mpsStrings = makeMPSstrings(mps,k-1); %needed for grad computation
        [optEnergy,gradEnergy] =   solveSDPandComputeGrad(k,W,h,rho_kp1,omega,D,d,mpsStrings,YalmipOps);
        grad= -gradEnergy; % minus sign because we minimize(-f) instead of max(f)
    else
        [optEnergy]      =   solveSDPandComputeGrad(k,W,h,rho_kp1,omega,D,d,[],YalmipOps);
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

function [optEnergy,grad,solDiagnost] = solveSDPandComputeGrad(k,W,h,rho_kp1,omega,D,d,mpsStrings,YalmipOps)

 
[constraints] = setConstraints(k,rho_kp1,omega,W,D,d);

 
objFunc=trace(kron(h,speye(d^(k-1))) * rho_kp1);
solDiagnost = optimize(constraints,objFunc,YalmipOps);

optEnergy = value( objFunc);

%% compute gradient w.r.t. mps tensor
if nargout >1

    gamma = struct('L',zeros(d*D),'R',zeros(d*D));
    
    fns = fieldnames(gamma);
    
    for fni=1:2
        yalmipVersion = yalmip('version');
        switch yalmipVersion
            case {'20210331'}
                G = dual(constraints(fns{fni})) ;
                gamma.(fns{fni}) = (G + G')/2   ;
            case {'20210609'}
                error('for some reason dual var is incorrectly read out in this version =?')
                % lowerTriangIndices = tril(true(d*D));
                % G = zeros(d*D);
                % G(lowerTriangIndices) = dual(constraints(fns{fni})) ;
                % gamma.(fns{fni}) = (G + G')/2   ;
        end
    end
    rhoOpt=value(rho_kp1);
    
    IxW=kron(speye(d),W);
    WxI=kron(W,speye(d));

    bondDim = sqrt(D);

    
    % TAKE d/dM OF trace( ...
    %    
    %    |    |             |
    %   GAMMA_L_GAMMA_L_GAMMA_L
    %    |    |             |
    %    |    .-M--M--M--M--.
    %    |      |  |  |  |
    %    RHORHORHORHORHORHO
    %    |      |  |  |  |
    %    |    .-M--M--M--M-.
    %    |    |            |
    % 
    % and the same with the mps on the left.
    %
    % we sum over all mps positions = 1..4
    % eg for pos=2
    % we need to "trace" the following (% are open inidices)
    %    |    |             |
    %   GAMMA_L_GAMMA_L_GAMMA_L
    %    |    |             |
    %    |    .-M--M--M--M--.
    %    |      |  |  |  |
    %    RHORHORHORHORHORHO
    %    |      |  |  |  |
    %    |    .-M-%%%-M--M-.
    %    |    |            |
    %    
    % we first contract the following
    %    1   -1            -2
    %    |    |             |
    %   GAMMA_L_GAMMA_L_GAMMA_L           -1     -2 
    %    |    |             |              |      |
    %    |    .-M--M--M--M--.       =  RsandwichRsandwich
    %    |      |  |  |  |                |  |  |  | 
    %    RHORHORHORHORHORHO              -3 -4 -5 -6
    %    |      |  |  |  |
    %    |     -3 -4 -5 -6
    %    1    
    %    
    Rsandwich = decomposeLegs(gamma.R * IxW * rhoOpt,[d, bondDim, bondDim],d*ones(1,k+1));
    Rsandwich = ncon({Rsandwich}, {[1, -1, -2, 1, -3:-1:-3-k+1]});
    
    Lsandwich = decomposeLegs(gamma.L * WxI * rhoOpt,[bondDim, bondDim, d],d*ones(1,k+1));
    Lsandwich = ncon({Lsandwich}, {[-1, -2, 1, -3:-1:-3-k+1, 1]});
    
    sandwich = Lsandwich + Rsandwich;  %from here on everything we do to L and R is the same so we can take the sum now.
    
    
    
    
    % RECALL: MPS INDEX CONVENTION
    %       -k-1          -k-2
    %        |             |
    %        .-M--M--M--M--.
    %          |  |  |  |
    %         -1 -2 .. -k
    
    % for position 1:
    %
    %                       -1      4=k
    %                        |      |
    %      (-2)             sandwichsandwich 
    %       |                |  |  |  | 
    % (-1)-dM-(-3) =        -2  1  2  3
    %                           1  2  3
    %                           |  |  | 
    %                      (-3)-M--M--M--4=k
    %                                  
    dM = ncon({sandwich, mpsStrings{k-1} },{[-1 k -2 1:k-1], [1:k-1 -3 k] });
    
    % for position pos=2:
    %
    %                        k      k+1
    %                        |      |
    %      (-2)             sandwichsandwich
    %       |                |  |  |  | 
    % (-1)-dM-(-3) =   pos-1=1 -2  2  3
    %                     1           2  3
    %                     |           |  | 
    %                   k-M-(-1) (-3)-M--M--k+1
    %                                  
    for pos=2:k-1
        dM = dM + ncon({sandwich,                          mpsStrings{pos-1},  mpsStrings{k-pos} },...
                       {[k, k+1, 1:pos-1, -2, pos:k-1], [1:pos-1 k -1], [pos:k-1 -3 k+1] },...
                       [k 1:pos-1 pos:k-1 k+1]); %contraction order
    end
    
    % for position k:
    %
    %                        k     (-3)
    %                        |      |
    %      (-2)             sandwichsandwich 
    %       |                |  |  |  | 
    % (-1)-dM-(-3) =         1  2  3 (-2)
    %                        1  2  3
    %                        |  |  | 
    %                      k-M--M--M--(-1)
    %                                  
    dM = dM + ncon({sandwich, mpsStrings{k-1} },{[k, -3, 1:k-1, -2], [1:k-1, k, -1] });
    
    grad= -2* dM(:); % minus sign because of SDP grad formula. factor 2 becuse of d(bra)/d... + d(ket)/d...
    
    
    
    
    % direct computation of gradient for debugging
    mps=permute(mpsStrings{1},[2 1 3]);
    
    indCell={[-(k+1), -1, 1]};
    for l=1:k-2
        indCell = horzcat(indCell, [l, -(l+1), l+1]);
    end
        indCell = horzcat(indCell, [(k-1), -(k) -(k+2)]); % add last entry

        
    dM2 = zeros(size(dM));
   
    % mps tensor with inds (i,j,l)
    for i=1:bondDim
        for j=1:d
            for l=1:bondDim
                Eijl  = zeros(bondDim,d,bondDim);
                Eijl(i,j,l) = 1;

                for pos=1:k
                    % multiply k mps together for Eijl in each position
                    toContract = [repmat({mps},1,pos-1), Eijl, repmat({mps},1,k-pos)] ;
                    Vijl_pos = reshape( ncon(toContract ,indCell,[],[-k-2, -k-1, -k:1:-1]),bondDim^2,d^k );


                    dM2(i,j,l) = dM2(i,j,l) + trace(...
                                                    gamma.L * kron(Vijl_pos,eye(d)) * rhoOpt * WxI' + ...
                                                    gamma.R * kron(eye(d),Vijl_pos) * rhoOpt * IxW' );
                end
            end
        end
    end
    grad2= -2*  dM2(:); % minus sign because of SDP grad formula. factor 2 becuse of d(bra)/d... + d(ket)/d...
    assert(norm(grad - grad2)<1e-14)
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



