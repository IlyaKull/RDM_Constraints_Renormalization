function [V0,L,R,L_kraus,R_kraus] = CoarseGrainingMaps(n,k0,vuMPS,CGops)
% construct coarse graining maps from input MPS for a chain of length n
% CGops specifies which maps to construct. The possible methods are:
    %       mpsLeft:
        %       useses the mps in LEFT gauge directly.
    %       isometriesLeft:
        %      	makes isometries from MPS in LEFT gauge. This can
        %       be done only up to a certain length because it involves a qr
        %       decomposition of an n-body map (from n spins to both bonds).  
        %       NOTE: the transformation of mps into isometries might involve an
        %       ill-conditioned map and in that case should be avoided as it is unstable
    %       mixLeft:
    %           start with isometries up to length defined by the variable CGops.isometriesUpTo_nmax 
    %           and from then on swich to mpsLeft. This specifying a correction
    %           matrix that relates the two methods.
    %       mpsFlip:
        %       useses the mps in the unconventional "FLIP" gauge, where the mps on the
        %       right is in LEFT gauge and the left one is in the RIGHT gauge. The
        %       reason for doing this is that this is the correct isometry direction
        %       which is needed when checking the feasibility of the dual point (see
        %       FeasibleDualPt.m)
        %       NOTE: the gauge transormation to this form might be ill-conditioned. In
        %       such case the use of this method is unstable 
    %       isometriesFlip:
    %           same as isometriesLeft only with Flip-gauged mps.
    %       mixFlip:
    %           same as mixLeft only with Flip-gauged mps.


if strncmp(CGops.whichMaps,'mix',3) 
    if  CGops.isometriesUpTo_nmax >= n
        nmax=n;
    else
        nmax =  CGops.isometriesUpTo_nmax ;
    end
    if nmax<=k0
       nmax=k0+1;
    end
end

fprintf('~~~~~~~ ')
switch  CGops.whichMaps
    case 'isometriesLeft'
        fprintf('Computing isometries left...\n')
        [V0,L,R]= isometriesFromMPS_QR(vuMPS.AL,k0,n);
    case 'isometriesFlip'
        fprintf('Computing isometries flip...\n')
        [V0,L,R]= isometriesFromMPS_QR_Flip(vuMPS,k0,n);
    case 'mpsLeft'
        fprintf('Computing MPS left...\n')
        [V0,L,R]= CGmapsFromMPS_Left(vuMPS.AL,k0,n);
    case 'mpsFlip'
        fprintf('Computing MPS flip...\n')
        [V0,L,R]= CGmapsFromMPS_Flip(vuMPS,k0,n);
    case 'mixLeft'  %from left mps
        fprintf('Mixed CG maps \n Computing isometries left...\n')
        [V0,Lqr,Rqr,lastP]= isometriesFromMPS_QR(vuMPS.AL,k0, nmax); %compute isometries up to nmax
        fprintf('Computing MPS left...\n')
        [~,Lmps,Rmps]= CGmapsFromMPS_Left(vuMPS.AL,k0,n); % in addition compute mps 
    case 'mixFlip'
        fprintf('Mixed CG maps \n Computing isometries flip...\n')
        [V0,Lqr,Rqr,lastP]= isometriesFromMPS_QR_Flip(vuMPS,k0, nmax); %compute isometries up to nmax
        fprintf('Computing MPS flip...\n')
        [~,Lmps,Rmps]= CGmapsFromMPS_Flip(vuMPS,k0,n); % in addition compute mps 
end
fprintf('~~~~~~~ Done computing CG maps. \n')

if strncmp(CGops.whichMaps,'mix',3) 
% the mixing of isometries and mps tensors for coarse graining means that
% we start with Lqr and Rqr isometries and at level k=nmax we switch to
% Lmps,Rmps CG maps. For this to be consistent we need to apply the 'lastP'
% matrix at level k=nmax-1
    if n==k0
        L=[];
        R=[];
    else
        L=cell(1,n-1);
        R=cell(1,n-1);
    end
    for k=k0:(nmax-2)
        L{k}=Lqr{k};
        R{k}=Rqr{k};
    end
    if n==nmax
        L{nmax-1}=Lqr{nmax-1}; %no need for the correction using lastP if this is the last level
        R{nmax-1}=Rqr{nmax-1};
    else   % n>nmax.. nmax cannot be bigger than n, see above
        L{nmax-1}=lastP*Lqr{nmax-1}; 
        R{nmax-1}=lastP*Rqr{nmax-1};
        for k = nmax:n-1
            L{k}=Lmps{k};
            R{k}=Rmps{k};
        end
    end
end
