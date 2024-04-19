function [V0 ,L,R,stats,isometries,MPSProdMatSave]= isometriesFromMPS_QR_DEBUG(vuMPS,k0,n,debug)
% construct isometries from uniform mps
if nargin==3
    debug=false;
end
    
% reshape vuMPS from cell to 3dim array 
d=length(vuMPS);
D=size(vuMPS{1,1},1);
mps=zeros(D,d,D);
for l=1:d
    mps(:,l,:)=vuMPS{l,1};
end
 
if debug
     MPSProdMatSave=cell(n,1);
     isometries=cell(n,1);
end
PmpsL=cell(n,1);
PmpsR=cell(n,1);
%% compute mps contractions
% this will not be needed for all k in the working version of the function
% (we only require the initial V and all Rs and Ls above it)
% here all Vs are computed in order to check the consictencey of the
% calculations 
  %for lower k0 P^2 has small negative eigenvals
for k=k0:n
    % compute the tensor composed of k mps tensors
    %contarction pattern: 
    %                    (-k-1) --A-(1)-A-(2)-A-(3)-A-(4)-A-- (-k-2)
    %                             |     |     |     |     |
    %                             -1    -2   -3    -4    -5
    %
    % because we want 
    %                   (-k-2)--VV--(-5)
    %                           VV--(-4)
    %                           VV--(-3)
    %                           VV--(-2)
    %                   (-k-1)--VV--(-1)
    
    indCell={[-(k+1), -1, 1]};
    for l=1:k-2
        indCell = horzcat(indCell, [l, -(l+1), l+1]);
    end
        indCell = horzcat(indCell, [(k-1), -(k) -(k+2)]); % add last entry
    
    MPSProd=ncon(repmat({mps},1,k),indCell,[],[-k-2, -k-1, -k:1:-1]);
    % reshape into a matrix from k physical spins into two bonds 
    MPSProdMat=reshape(MPSProd,D^2,d^k); 
    %to obtain an isometry from  MPSProdMat   P^-1 needs to be applied to it
    % where P^2 is the following D^2xD^2 matrix:
    [Q,r]=qr(MPSProdMat',0);
    P{k}=r';
    if k==k0
        V0=Q';
        L=[];
        R=[];
    end
    if debug  %save isometries only in debug mode
        MPSProdMatSave{k}=MPSProdMat;
        isometries{k}=Q';
    end
  
end
%% test isometries  
if debug
    stats.normVV_I=NaN(n,1);
    for k=k0:n
     stats.normVV_I(k)= norm(isometries{k}*isometries{k}'-eye(D^2)) ;
    end
end
 %% Left and Right Isometries
%  the left and right isometries should satisfy:
%    L_k *(I,V_k)  = V_k+1  ;    R_k *(V_k,I)  = V_k+1 
%
%   --LL--VV--   --VV--     ;   --RR------   --VV-- 
%     LL  VV--     VV--     ;     RR--VV--     VV-- 
%     LL  VV--  =  VV--     ;     RR  VV--  =  VV-- 
%     LL--VV--     VV--     ;     RR  VV--     VV-- 
%   --LL------   --VV--     ;   --RR--VV--   --VV-- 

% i.e.
%     L_k = (P_k+1)^-1 * (PmpsL)_k
%
%   --LL--     --Pinv-------PP--    %   --RR----    --Pinv   .-------
%     LL         Pinv       PP      %     RR--        Pinv   |   PP--      
%     LL     =   Pinv--mmm--PP      %     RR     =    Pinv--www--PP      
%     LL--       Pinv   |   PP--    %     RR          Pinv       PP      
%   --LL----   --Pinv   .-------    %   --RR--      --Pinv-------PP--     
%
% where  P is r' from the qr decomposition and we define (PmpsL)_k and (PmpsR)_k :
%
%    -------PP--      --PmpsL--    ;     .-------        PmpsL---
%           PP          PmpsL      ;     |   PP--        PmpsL--
%    --mps--PP     =: --PmpsL      ;  --www--PP     =: --PmpsL
%       |   PP--        PmpsL--    ;         PP          PmpsL
%       .-------        PmpsL---   ;  -------PP--      --PmpsL--
% 
% (mmm) and (www) stand for the mps tensor  :
%
%                             (d) 
%   (l)--mmm--(r)              |  
%         |         =    (r)--www--(l)
%        (d)                
% as we don't want to compute r^-1 explicitly
% we solve  
%   (r_k+1) * L_k = (PmpsL)_k 
% for L_k, and similarly for R_k:
%   (r_k+1) * R_k = (PmpsR)_k


% first compute PmpsL and PmpsR
%

mpsInflatedR=zeros(D^2,d,D^2);
mpsInflatedL=zeros(D^2,d,D^2);
for l=1:d
    mpsInflatedR(:,l,:)=tensor(eye(D),vuMPS{l,1});
    mpsInflatedL(:,l,:)=tensor(vuMPS{l,1},eye(D));
end
%  PmpsL:
% 
%  (-1)------(1)-PP--(-3)     
%  (-1)--mps-(1)-PP--(-3)
%         |           
%         .----------(-2)       
% 
%  PmpsR:
%         .----------(-3)
%         |
%  (-1)--www-(1)-PP--(-2)     
%  (-1)------(1)-PP--(-2)           
%          
for k=k0:n     
    % compute PL{k} and PR{k}
   PmpsL{k}=ncon({mpsInflatedL,P{k}},{[-1,-3,1],[1,-2]})  ;
   PmpsL{k}=reshape(PmpsL{k},D^2,(D^2)*d); %reshape to a matrix
   PmpsR{k}=ncon({mpsInflatedR,P{k}},{[1,-2,-1],[1,-3]}) ; 
   PmpsR{k}=reshape(PmpsR{k},D^2,(D^2)*d); %reshape to a matrix
 
end

%% test PmpsL & PmpsR
if debug
    stats.maxPmpsL_MPSprod=NaN(n-1,1);
    stats.maxPmpsR_MPSprod=NaN(n-1,1);
    for k=k0:n-1   
        stats.maxPmpsL_MPSprod(k)=max(abs(...
               PmpsL{k}*tensor(eye(2),isometries{k}) - MPSProdMatSave{k+1}...
               ),[],'all');
        stats.maxPmpsR_MPSprod(k)=max(abs(...
               PmpsR{k}*tensor(isometries{k},eye(2)) - MPSProdMatSave{k+1}...
               ),[],'all');
    end
     stats.maxPmpsL_PmpsR=NaN(n,1);
    for k=k0:n
        stats.maxPmpsL_PmpsR(k)= max(abs(PmpsL{k}*PmpsL{k}' - PmpsR{k}*PmpsR{k}'),[],'all');
    end
end
%% compute left and right isometries
% recall that L and R are defined by P{k+1}*L{k} = MPS*r{k} =:PmpsL{k}
% ie  P{k+1}*L{k} = PL{k}
% and P{k+1}*R{k} = PR{k}
% we solve these equations for L and R
%  linsolve: X = linsolve(A,B) solves the linear system AX = B
recipCondOut=nan(n,1);
for k=k0:n-1
        linsolveOps.LT=true; %lower triangular.  
    [Ltemp,recipCondOut(k+1)]=linsolve(P{k+1} ,PmpsL{k} ,linsolveOps);
    L{k}=Ltemp;
    Rtemp=linsolve(P{k+1} ,PmpsR{k} ,linsolveOps);
    R{k}=Rtemp;
end

stats.conditionNumbersOfP=recipCondOut.^-1;
%%   test 
if debug
    stats.normLL_I=NaN(n-1,1);
    stats.normRR_I=NaN(n-1,1);
    stats.normRV_LV=NaN(n-1,1);
    stats.normRV_V=NaN(n-1,1);
    stats.normLV_V=NaN(n-1,1);
    for k=k0:n-1
        stats.normLL_I(k)=    norm(L{k}*L{k}'-eye(D^2));
        stats.normRR_I(k)=    norm(R{k}*R{k}'-eye(D^2));
        stats.normRV_LV(k) =  norm(R{k}*tensor(isometries{k},eye(d))-...
            L{k}*tensor(eye(d),isometries{k}) );
        stats.normRV_V(k) =  norm(isometries{k+1}-...
            R{k}*tensor(isometries{k},eye(d)));
        stats.normLV_V(k) =  norm(isometries{k+1}-...
            L{k}*tensor(eye(d),isometries{k})) ;
    end 
end
