function [] = warmStartSCSMF(fileName,maxIter,typ,lookback,safeguard)
    

    dataFSpath=getenv('DATA'); %save files to DATA file system when working on cluster
    if isempty(dataFSpath)
        savePathAndFilename = fileName;
    else
        savePathAndFilename = [dataFSpath,'/translInvSpinChain/outFiles/matFiles/',fileName];
    end
    
    load(savePathAndFilename,'PD')
    PD.maxIter = maxIter;
    
    if ~isfield(PD,'model')
        modelNames={'heisenberg2_','heisenberg2U_','heisenberg3_','XXZ_','XXZU_','XY_'};
        ind=zeros(length(modelNames),1);
        for j=1:length(modelNames) 
           ind(j)= contains(fileName,modelNames{j});
        end
        PD.model=modelNames{find(ind)}(1:end-1);
    end
    
    PD.AAtyp=typ;    
    PD.AAlookback=lookback; 
    PD.AAsafeguard=safeguard; 
    
    [energy,~,~,~,~,~,PD,AAmisses]  = runSCS_AA(PD.model,PD.D,PD.n,PD);
    
    format long
    fprintf('===========================================================================================================\n')
    fprintf('===========================================================================================================\n')
    fprintf('===========================================================================================================\n')
     fprintf('AA type: %d , lookback: %d , safegaurd: %d \n', typ,lookback,safeguard)
    fprintf('Acceleration misses: % d (out of %d total AA steps) \n',AAmisses, PD.maxIter/lookback)
    disp('Lower bound energy difference history:')
    disp(array2table([energy(:,1) , (PD.Eexact - energy(:,2))],'VariableNames',{'iteration','Delta energy'}))
end