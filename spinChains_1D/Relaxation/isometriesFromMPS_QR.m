function [V0 ,L,R,lastP,condNumsOut]= isometriesFromMPS_QR(vuMPS,k0,n)
% construct isometries from uniform mps
% this function is not commented at all! for detailed documentation see:
% isometriesFromMPS_QR_DEBUG.m
condNumsOut=[];
  
    
% reshape vuMPS from cell to 3dim array 
d=length(vuMPS);
D=size(vuMPS{1,1},1);
mps=zeros(D,d,D);
for l=1:d
    mps(:,l,:)=vuMPS{l,1};
end
 
% create 'inflated mps' for later use
mpsInflatedR=zeros(D^2,d,D^2);
mpsInflatedL=zeros(D^2,d,D^2);
for l=1:d
    mpsInflatedR(:,l,:)=tensor(eye(D),vuMPS{l,1});
    mpsInflatedL(:,l,:)=tensor(vuMPS{l,1},eye(D));
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
indCell={[-(k0+1), -1, 1]};
for l=1:k0-2
    indCell = horzcat(indCell, [l, -(l+1), l+1]);
end
    indCell = horzcat(indCell, [(k0-1), -(k0) -(k0+2)]); % add last entry

MPSProd=ncon(repmat({mps},1,k0),indCell,[],[-k0-2, -k0-1, -k0:1:-1]);
% reshape into a matrix from k physical spins into two bonds 
MPSProd=reshape(MPSProd,D^2,d^k0); 

% first isometry and P matrix:
[Q,r]=qr(MPSProd',0);
V0=Q';
clear Q
P=r';
clear r
if k0==n
    return
end

PmpsL=ncon({mpsInflatedL,P},{[-1,-3,1],[1,-2]})  ;
PmpsL=reshape(PmpsL,D^2,(D^2)*d); %reshape to a matrix
PmpsR=ncon({mpsInflatedR,P},{[1,-2,-1],[1,-3]}) ; 
PmpsR=reshape(PmpsR,D^2,(D^2)*d); %reshape to a matrix

 
for k=(k0+1):n
    
    % add one more mps tensor to MPS product
    MPSProd=ncon({mpsInflatedL,MPSProd},{[-1,-3,1],[1,-2]})  ;
    MPSProd=reshape(MPSProd,D^2,d^k); %reshape to a matrix
        
    % obtain only triangular r from qr
    [~,r]=qr(MPSProd',0);
    P=r';
    % L  and R are computed from current P and previous PmpsL/R. see
    % documentation in isometriesFromMPS_QR_DEBUG
    condP=cond(P);
    fprintf('condition number of P%d: %0.3g \n',k, condP)
    condNumsOut=[condNumsOut,condP];
    
    linsolveOps.LT=true; %options for linear solver: LT=lower triangular.
    L{k-1}=linsolve(P,PmpsL ,linsolveOps);
    R{k-1}=linsolve(P,PmpsR,linsolveOps);
    
    PmpsL=ncon({mpsInflatedL,P},{[-1,-3,1],[1,-2]})  ;
    PmpsL=reshape(PmpsL,D^2,(D^2)*d); %reshape to a matrix
    PmpsR=ncon({mpsInflatedR,P},{[1,-2,-1],[1,-3]}) ; 
    PmpsR=reshape(PmpsR,D^2,(D^2)*d); %reshape to a matrix
  
end
 lastP=P;
