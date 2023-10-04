clearvars
sss= tic 
 

modelStr = 'heisenberg2U'
n = 12;
D = 5;
f1=figure;
e1=axes(f1);
f2=figure;
e2=axes(f2);

ops.maxIter=1000;
ops.iterStride=100;
ops.saveProgress=true;
ops.saveStride=500;
ops.parallelize='single';

[energy,~,~,DeltaLower,~,res,PD,misses]  = runSCS_AA(modelStr,D,n,ops);


%%
legendCell={};

% % [energyV,~,~,DeltaLowerV,~,resV,PD]  = runSCS(modelStr,D,n,ops);
% plot( e1 ,energyV(:,1),vecnorm(resV,2,2))
% hold(e1,'on')
% plot( e2,energyV(:,1),ones(size(energyV(:,1) ))*DeltaLowerV,'--')
% hold(e2,'on')
% legendCell={'vanilla'};

for typ=[ 1  ]
    for lookback=[ 10  ]
        for safeguard=[ 3 ]

            ops.AAtyp=typ;    
            ops.AAlookback=lookback; 
            ops.AAsafeguard=safeguard; 

            [energy,~,~,DeltaLower,~,res,PD,misses]  = runSCS_AA(modelStr,D,n,ops);
            fprintf('AA type: %d , lookback: %d , safegaurd: %d \n', typ,lookback,safeguard)
            fprintf('Acceleration misses: % d (out of %d total AA steps) \n',misses, ops.maxIter/lookback)
            disp(array2table([energy(:,1) , (PD.Eexact - energy(:,2))],'VariableNames',{'iteration','DeltaE'}))

            legendCell=[legendCell,['typ=',num2str(typ),' mem=',num2str(lookback),' safeguard=',num2str(safeguard),' miss ratio=',num2str(misses/ops.maxIter*lookback)]];

            plot(e1, energy(:,1),vecnorm(res,2,2))
            plot( e2,energy(:,1),ones(size(energy(:,1) ))*DeltaLower,'--')
            
%              
        end
    end
end
toc(sss)
legend(e1,legendCell)
legend(e2,legendCell)

 set(e1,'XScale','linear','YScale','log')
 set(e2,'XScale','linear','YScale','log')
% norm(energy1-energy2)
% norm(reportTimes1 - reportTimes2)
% norm(DeltaLower1 - DeltaLower2)
% norm(res1-res2)