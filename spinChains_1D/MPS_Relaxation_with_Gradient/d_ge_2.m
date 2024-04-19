
d=2;
rng(1,'twister')
h = randn(d^2); h=h+h';

H.interactionTerm = h;
H.localDim = d;

for N=3:6
    
    D=d^(N-2);
    [energyFinalLTI(N),energyInitLTI(N)]=LTI_grad_opt_test(D,N,H);
end


%%
filename='one_W_d_eq2_rnadnH.mat';

DeltaLTI_final=cell(6,1);
DeltaLTI_init=cell(6,1);

for N=5:6
    for seed=1:100
        for D = 2:5
            [energyFinal(D,seed),energyInit(D,seed),Wopt{D,seed},W0{D,seed},exitflag{D,seed},output{D,seed}]=...
                LTI_grad_opt_test(D,N,H,seed);
        end
    end
     


DeltaLTI_final{N}=energyInitLTI(N) - energyFinal;
DeltaLTI_init{N}=energyInitLTI(N) - energyInit;


figure
markers='ox^sp';
for D=2:5
    scatter(DeltaLTI_init{N}(D,:),DeltaLTI_final{N}(D,:),markers(D))
    hold on
end
ax=gca;
ax.YScale='log';
ax.XScale='log';
title(['N=',num2str(N)])
  save(filename)

end
