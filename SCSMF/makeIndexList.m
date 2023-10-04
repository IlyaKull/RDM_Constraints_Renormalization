function [xinds,yinds]=makeIndexList(d,D,k0,n)
% xinds is the index list for the dual variables 
% for example
% %    
%                 indsFromTo     dims
%                ____________    ____
% 
%     epsil         1       1      1 
%     alpha         2      65      8 
%     betaL        66     389     18 
%     betaR       390     713     18 
%     gammaL3     714    1037     18 
%     gammaR3    1038    1361     18 
%     gammaL4    1362    1685     18 
%     gammaR4    1686    2009     18 
%     gammaL5    2010    2333     18 
%     gammaR5    2334    2657     18 
%
% yinds is the index list for the dual variables 
% for example

%% xinds
varNames={'epsil';'alpha';'betaL';'betaR'};

dimAlpha=d^k0;
dimBeta=d*D*D;
dims=[1 ; dimAlpha; dimBeta; dimBeta];

for i=k0:n-1
    varNames=[varNames;{['gammaL',num2str(i)];['gammaR',num2str(i)]}];
    dims=[dims;dimBeta;dimBeta];
end

indsFromTo=nan(length(varNames),2);
indsFromTo(1,1)=1; indsFromTo(1,2)=1;               %indices of epsil var
indsFromTo(2,1)=2; indsFromTo(2,2)=2+dimAlpha^2-1; %indices of alpha var
for i=3:length(varNames)
     indsFromTo(i,1)=indsFromTo(i-1,2)+1;
     indsFromTo(i,2)=indsFromTo(i,1)+dimBeta^2-1;
end
varNum=[1:length(varNames)]';
xinds=table(varNum,indsFromTo,dims,'RowNames',varNames);
% to acess table entries use xinds{'epsil',:} 
% or xinds('epsil',:).Variables

%% yinds
varNames={'rho'};

dimRho=d^(k0+1);
dimSigma=d*D*D*d;
dims=dimRho ;

for i=k0:n
    varNames=[varNames;['sigma',num2str(i)]];
    dims=[dims;dimSigma];
end

indsFromTo=nan(length(varNames),2);
indsFromTo(1,1)=1; indsFromTo(1,2)=dimRho^2;               %indices of rho var
for i=2:length(varNames)
     indsFromTo(i,1)=indsFromTo(i-1,2)+1;
     indsFromTo(i,2)=indsFromTo(i,1)+dimSigma^2-1;
end
varNum=[1:length(varNames)]';
yinds=table(varNum,indsFromTo,dims,'RowNames',varNames);


