function [V0 ,L,R,mps]= CGmapsFromMPS_Left(vuMPS,k0,n,inputFormat)
% instead of computing isometries from MPS JUST USE MPS products
% this makes an SDP which is equivalent to the one that uses isometries
% when the MPS is injective. This is because in this case the isometries
% and are related to the MPS products by an invertible map (P in the
% previous functions)
    

if nargin<4
    inputFormat=0;
end

if inputFormat==0
    % reshape vuMPS from cell to 3dim array 
    d=length(vuMPS);
    D=size(vuMPS{1,1},1);
    mps=zeros(D,d,D);
    for l=1:d
        mps(:,l,:)=vuMPS{l,1};
    end

else
   mps=vuMPS;
   D=size(mps,1);
   d=size(mps,2);
   assert(D==size(mps,3),'CHECK MPS DIMS')
end

% create 'inflated mps' for later use
mpsInflatedR=zeros(D^2,d,D^2);
mpsInflatedL=zeros(D^2,d,D^2);
for l=1:d
    mpsInflatedR(:,l,:)=tensor(eye(D),squeeze(mps(:,l,:)));
    mpsInflatedL(:,l,:)=tensor(squeeze(mps(:,l,:)),eye(D));
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
V0=reshape(MPSProd,D^2,d^k0); 
 
%%

if k0==n
    return
end
 
for k=k0:n-1
    %ncon is just used for reshaping into the corect order here 
    L{k}=reshape(ncon({mpsInflatedL},{[-1,-3,-2]}) ,D^2,(D^2)*d);  
    R{k}=reshape(ncon({mpsInflatedR},{[-3,-2,-1]}) ,D^2,(D^2)*d); 
end
 