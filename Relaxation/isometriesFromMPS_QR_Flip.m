function [V0 ,L,R,lastP,condNumsOut]= isometriesFromMPS_QR_Flip(vuMPS,k0,n)
% construct isometries from uniform mps
% this function is not commented at all! for detailed documentation see:
% isometriesFromMPS_QR_DEBUG.m
condNumsOut=[];
linsolveOps.LT=true; %options for linear solver: LT=lower triangular.  
   
% the 'Flip' version of the function makes coarse graining maps
% -AR-AR-...-Cinv-AL-...-AL-
% ie the left tensor is on the right and the right one is on the left


% reshape vuMPS from cell to 3dim array 
d=length(vuMPS.AL);
D=size(vuMPS.AL{1,1},1);
mpsL=zeros(D,d,D);
mpsCinv=zeros(D,d,D);
mpsR=zeros(D,d,D);
testAR = nan(d,1);
testAL = nan(d,1);
testAC = nan(d,1);

for jj=1:d
    mpsL(:,jj,:)=vuMPS.AL{jj,1};
    mpsR(:,jj,:)=vuMPS.AR{1,jj};
    % solve for AL and AR simultaneusly
    [mpsCinv(:,jj,:)]= solveForLandR(vuMPS.C, vuMPS.AL{jj,1} ,vuMPS.AR{1,jj}) ;
    %tests
    testAR(jj)=max(max(abs(squeeze(mpsCinv(:,jj,:))*vuMPS.C - vuMPS.AR{jj})));
    testAL(jj)=max(max(abs(   vuMPS.C*squeeze(mpsCinv(:,jj,:)) - vuMPS.AL{jj})));
    testAC(jj)=max(max(abs(vuMPS.C*squeeze(mpsCinv(:,jj,:))*vuMPS.C - vuMPS.AC{jj})));
 end
fprintf('condition number of vuMPS.C: %0.3g \n', cond(vuMPS.C))
fprintf('|ACinv*C-AR|_max= %0.3g \n', max(testAR))
fprintf('|C*ACinv-AL|_max= %0.3g \n', max(testAL))
fprintf('|C*ACinv*C-AC|_max= %0.3g \n', max(testAC))

% create 'inflated mps' (tensor an identity on the bond)
% this will be reshaped in the end to form the CG map
mpsInflatedR=zeros(D^2,d,D^2);
mpsInflatedL=zeros(D^2,d,D^2);
for l=1:d
    mpsInflatedR(:,l,:)=tensor(eye(D),vuMPS.AL{l,1}); %!!! in the 'flip' version AR is on the Left and viceverca
    mpsInflatedL(:,l,:)=tensor(vuMPS.AR{1,l},eye(D));
end

% no left and right isometries if there is only one coarser-graining step
if n==k0
    L=[];
    R=[];
else
    L=cell(1,n-1);
    R=cell(1,n-1);
end
 

% multiply k0 mps together
% the product is always of the form -AR-AR-...-ACinv-AL-...-AL-
    % indices for ncon:
if k0==1
    indCell= {[-2,-1,-3]};
else
    indCell={[-(k0+1), -1, 1]};
    for l=1:k0-2
        indCell = horzcat(indCell, [l, -(l+1), l+1]);
    end
        indCell = horzcat(indCell, [(k0-1), -(k0) -(k0+2)]); % add last entry
end

if k0==1
    mpsToMultiply=mpsCinv;
elseif mod(k0,2)==0
    mpsToMultiply= horzcat(repmat({mpsR},1,k0/2),mpsCinv,repmat({mpsL},1,(k0/2)-1));
else 
    mpsToMultiply= horzcat(repmat({mpsR},1,(k0-1)/2),mpsCinv,repmat({mpsL},1,(k0-1)/2));
end

MPSProd=ncon(mpsToMultiply,indCell,[],[-k0-2, -k0-1, -k0:1:-1]);
% reshape into a matrix from k physical spins into two bonds 
MPSProd=reshape(MPSProd,D^2,d^k0); 

% first isometry and P matrix:
[Q,r]=qr(MPSProd',0);
V0=Q'; clear Q
P=r'; clear r
if k0==n
    return
end



PmpsL=ncon({mpsInflatedL,P},{[-1,-3,1],[1,-2]})  ;
PmpsL=reshape(PmpsL,D^2,(D^2)*d); %reshape to a matrix
PmpsR=ncon({mpsInflatedR,P},{[1,-2,-1],[1,-3]}) ; 
PmpsR=reshape(PmpsR,D^2,(D^2)*d); %reshape to a matrix

 
for k=(k0+1):n
    
    % add one more mps tensor to MPS product alternating between left and
    % right
    if mod(k,2)
        MPSProd=ncon({mpsInflatedL,MPSProd},{[-1,-3,1],[1,-2]})  ;
    else
        MPSProd=ncon({mpsInflatedR,MPSProd},{[1,-2,-1],[1,-3]})  ;
    end
    MPSProd=reshape(MPSProd,D^2,d^k); %reshape to a matrix
        
    % obtain only triangular r from qr
    [~,r]=qr(MPSProd',0);
    P=r';
    % L  and R are computed from current P and previous PmpsL/R. see
    % documentation in ..._DEBUG
    condP=cond(P);
    fprintf('condition number of P%d: %0.3g \n',k, condP)
    condNumsOut=[condNumsOut,condP];
    L{k-1}=linsolve(P,PmpsL ,linsolveOps);
    R{k-1}=linsolve(P,PmpsR,linsolveOps);
    
    PmpsL=ncon({mpsInflatedL,P},{[-1,-3,1],[1,-2]})  ;
    PmpsL=reshape(PmpsL,D^2,(D^2)*d); %reshape to a matrix
    PmpsR=ncon({mpsInflatedR,P},{[1,-2,-1],[1,-3]}) ; 
    PmpsR=reshape(PmpsR,D^2,(D^2)*d); %reshape to a matrix
  
end
 lastP=P;
%% fix isometries (set singular values to exactly 1)
% for j=k0:n-1
%     [U,~,W]=svd(L{j});
%     L{j} = U*([eye(D^2),zeros(D^2,(d-1)*D^2)])*W';
%     [U,~,W]=svd(R{j});
%     R{j} = U*([eye(D^2),zeros(D^2,(d-1)*D^2)])*W';
% end