clear


for D=2:4 
    for N= 12:8:76
        t=tic;
        try
            [energyFinal(N,D),energyInit(N,D),mpsOpt{N,D},mpsInit{N,D},optOutput{N,D}]=MPSgrad_opt_full(D,N,0.1 ); 
        catch ms
            ERRS{N,D}  = ms;
        end    
        DeltaInit(N,D) = 0.25-log(2) -energyInit(N,D), 
        DeltaFinal(N,D) = 0.25-log(2) -energyFinal(N,D),
        totalTime(N,D) = toc(t);
        save('mps_grad_run.mat')
    end
end
%%
for D=5:7
    for N=60:20:100
        t=tic;
        try
            [energyFinal(N,D),energyInit(N,D),mpsOpt{N,D},mpsInit{N,D},optOutput{N,D}]=MPSgrad_opt_full(D,N,0.1 ); 
        catch ms
            ERRS{N,D}  = ms;
        end    
        DeltaInit(N,D) = 0.25-log(2) -energyInit(N,D), 
        DeltaFinal(N,D) = 0.25-log(2) -energyFinal(N,D),
        totalTime(N,D) = toc(t);
        save('mps_grad_run.mat')
    end
end
