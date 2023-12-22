function [V0 ,L,R,L_kraus,R_kraus]= CGmapsFromMPS_flip(vuMPS,k0,n)
% instead of computing isometries from MPS JUST USE MPS products
% this makes an SDP which is equivalent to the one that uses isometries
% when the MPS is injective. This is because in this case the isometries
% and are related to the MPS products by an invertible map (P in the
% previous functions)
    
% the 'Flip' version of the function makes coarse graining maps
% -AR-AR-...-Cinv-AL-...-AL-
% ie the left tensor is on the right and the right one is on the left
% the MPScenter map is related to this one by CxC on the bonds:
% -L-L-L-C-R-R-R-=-C-R-R-R-Cinv-L-L-L-C-    (LC=CR => L=C-R-Cinv,R=Cinv-L-C)


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
    %     mpsCinv(:,jj,:)= linsolve(vuMPS.C,vuMPS.AL{jj,1}); %%%%%%%%%%%%% ACinv=Cinv AC Cinv= Cinv AL= AR Cinv
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
for jj=1:d
    mpsInflatedR(:,jj,:)=tensor(eye(D),vuMPS.AL{jj,1}); %!!! in the 'flip' version AR is on the Left and viceverca
    mpsInflatedL(:,jj,:)=tensor(vuMPS.AR{1,jj},eye(D));
end

% no left and right isometries if there is only one coarser-graining step
if n==k0
    L=[];
    R=[];
     L_kraus=[];
    R_kraus=[];
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
    for jj=1:k0-2
        indCell = horzcat(indCell, [jj, -(jj+1), jj+1]);
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
V0=reshape(MPSProd,D^2,d^k0); 

if k0==n
    return
end
 
for k=k0:n-1
    %ncon is just used for reshaping into the corect order here 
    L{k}=reshape(ncon({mpsInflatedL},{[-1,-3,-2]}) ,D^2,(D^2)*d);  
    R{k}=reshape(ncon({mpsInflatedR},{[-3,-2,-1]}) ,D^2,(D^2)*d);  
end
  

%output the mps as kraus operators of a channel from (dxD) to (D)
% right and left are interchanged with the 'flip' options
L_kraus = reshape(permute(mpsR,[1,3,2]),D,d*D);
R_kraus = reshape(permute(mpsL,[3 2 1]),D,d*D);


