clear
d=2;
energyFinalLTI = nan(9,1);
energyInitLTI= nan(9,1);


d=2;
h = GetTwoSiteH([-1,-1,-Jz,0,0],d);     

H.interactionTerm = h;
H.localDim = d;

for N=4:7
    
    D=d^(N-2);
    [energyFinalLTI(N),energyInitLTI(N)]=LTI_grad_opt_test(D,N,H);
end
Wopt= cell(9,6,100);
W0= cell(9,6,100);
%%

for N=5 
    for seed=1:1000
        for D = 1 
            [energyFinal(N,D,seed),energyInit(N,D,seed),Wopt{N,D,seed},W0{N,D,seed},exitflag{N,D,seed},output{N,D,seed}]=...
                LTI_grad_opt_test(D,N,H,seed);
        
        end
    end
    save('one_W_random_seed.mat')
end
 

%%
energyFinal5=squeeze(energyFinal(5,:,:));
energyInit5=squeeze(energyInit(5,:,:));
DeltaLTI_final_5=energyInitLTI(5) - energyFinal5;
DeltaLTI_init_5=energyInitLTI(5) - energyInit5;
Wopt5=squeeze(Wopt(5,:,:));
W05=squeeze(W0(5,:,:));
exitflag5=squeeze(exitflag(5,:,:));
output5=squeeze(output(5,:,:));

%%
figure
markers='ox^sp';
for D=2:5
    scatter(DeltaLTI_init_5(D,:),DeltaLTI_final_5(D,:),markers(D))
    hold on
end
ax=gca;
ax.YScale='log';
ax.XScale='log';

