function [PD]=setProblemData(modelStr,D,n,ops)
 
modelInfo.modelName=modelStr;


%% method
% this speifies which coarse-graining maps to build from the MPS
% depending on the model, the MPS can be ill conditioned and therefore some
% methods do not perform as well as others
% 'mpsLeft' is always the safe choice

if isfield(ops,'method')
    methodOps.whichMaps=ops.method;
else
    switch modelInfo.modelName
            case {'TFI','heisenberg3','heisenberg2U','XXZU'}
                methodOps.whichMaps='mixFlip';
            otherwise
                methodOps.whichMaps='mpsLeft';
    end
end
methodOps.isometriesUpTo_nmax=12;

 

%% load saved MPS
% you can specify a file where the MPS is saved
% I had some default files specified 

if isfield(ops,'saved_MPS_filename')
    saved_MPS_filename = ops.saved_MPS_filename ;            
else
    if strcmp(modelInfo.modelName,'TFI') %g=1!
        saved_MPS_filename = 'savedMPScriticalTFI.mat';
    else
        if D>8
            saved_MPS_filename = 'savedMPStensorsD9t12.mat';
        else
            saved_MPS_filename = 'savedMPSwithSubLatRot.mat';
        end
    end
end


load(saved_MPS_filename,'Dmax','models','H_save','Eexact_save','vuMPS_save',...
            'Eupper_save','converged','xi_save');
 % the following just loads the variables needed from the saved file
 
        switch modelInfo.modelName
            case 'TFI'
%                   modelNameInds=find(cellfun(@(ss) strcmp(ss,modelInfo.modelName),models(:,1)));
%                 gInds=find(modelInfo.g == [models{:,2}]);
%                 indToLoad=intersect(gInds,modelNameInds);
                    indToLoad=find(cellfun(@(ss) strcmp(ss,modelInfo.modelName),models(:,1)));
            otherwise
                indToLoad=find(cellfun(@(ss) strcmp(ss,modelInfo.modelName),models(:,1)));
        end
        
        if D>Dmax || isempty(indToLoad) || ~converged(indToLoad,D)
            error('Cannnot find given model and D in file or VUMPS did not convege. No saved MPS.')
        else
            H= H_save{indToLoad,D};
            Eexact= Eexact_save(indToLoad,D);
            if strncmp(modelStr,'XXZ',3)
                Eexact = -0.617222045975845; %calculated with VUMPS with D=200
            end
            vuMPS= vuMPS_save{indToLoad,D};
            Eupper= Eupper_save(indToLoad,D);
            xi= xi_save(indToLoad,D);
            VUMPSconverged=converged(indToLoad,D); %in some models vumps doesn't converger for all Ds 
        end
        
d=sqrt(size(H,1));
k0=find(d.^(1:n)>D^2,1); % find the size of the first state 
if isempty(k0)
    k0=n;
end
   

% make the coarse graining maps from the MPS
[V0,L,R] = CoarseGrainingMaps(n,k0,vuMPS,methodOps);
 
%% make index tables 
% the optimization variables are vectors consisting of all the states
% \rho_0, \omega_1, \omega_2... 
% in order to perform operations on the states, their position within the
% large vectors have to be specified. Those indices are saved in the
% following table
% 
[XXinds,YYinds]=makeIndexList(d,D,k0,n);
xinds=XXinds{:,1:3}; %convert from table to matrix
yinds=YYinds{:,1:3};


% PD is the 'problem data' structure. it contains all the variables that
% need to be passed between functions 
%(yes, I also thing that an object oriented approach would have been nicer...)

PD=struct(      'model',modelStr,...
                'xinds',xinds,...
                'yinds',yinds,...
                'd',d,...
                'D',D,...
                'k0',k0,...
                'n',n,...
                'H',H,...
                'methodOps',methodOps ...
                );
                

%% SCS options 

% conjugate gradient initial precision
if isfield(ops,'CGtol')
    PD.CGtol = ops.CGtol;            
else
    PD.CGtol = 1e-4;            
end

% minimum precision (will not decrease below this value)
if isfield(ops,'CGtolMin')
    PD.CGtolMin = ops.CGtolMin;            
else
    PD.CGtolMin = 1e-12;            
end

% this setting changes the CG tolerance as a function of the current
% itration 
if isfield(ops,'CGtolFunc')
    PD.CGtolFunc = ops.CGtolFunc;            
else
    PD.CGtolFunc = @(x) max(PD.CGtol * x.^(-1.1) , PD.CGtolMin);            
end

%conj grad max iters
if isfield(ops,'CGmaxIter')
    PD.CGmaxIter = ops.CGmaxIter;            
else
    PD.CGmaxIter = 1e4;            
end

% conjugate gradient target precision for initial solve of Minvh (this is
% done once and the result is reused at every iteration)
if isfield(ops,'MinvhTol')
    PD.MinvhTol = ops.MinvhTol;            
else
    PD.MinvhTol = 1e-15;            
end

%conj grad max iters for initial solve of Minvh
if isfield(ops,'MinvhMaxIter')
    PD.MinvhMaxIter = ops.MinvhMaxIter;            
else
    PD.MinvhMaxIter = 1e5;            
end

% over relaxation parameter (internal SCS parameter, see paper)
if isfield(ops,'q')
    PD.q = ops.q;            
else
    PD.q = 1.5;            
end

% target tolerance for SCS residuals: [primalRes, dualRes, gapRes]
if isfield(ops,'residTol')
    PD.residTol = ops.residTol;            
else
    PD.residTol = [1 1 1]*1e-8;       
end

% SCS max iterations
if isfield(ops,'maxIter')
    PD.maxIter = ops.maxIter;            
else
    PD.maxIter = 2000;       
end


 % SCS iteration stride (check if converged every # iterations)
if isfield(ops,'iterStride')
    PD.iterStride = ops.iterStride;            
else
    PD.iterStride = 50;       
end

% % warning for deviation of matrices from being symmetric (during project to PSD
% % cone step)
% if isfield(ops,'asymmThresh')
%     PD.asymmThresh = ops.asymmThresh;            
% else
%     PD.asymmThresh = 1e-12;       
% end




% rescaling of c and b (see https://web.stanford.edu/~boyd/papers/pdf/scs_long.pdf (5))
% sigma and rho are  internal SCS parameters 
% sigma multiplies the Hamiltonian 
% choosing sigma=1/Normb makes the new b (b^hat) have norm 1
Normb=norm(H(:)) *2^( (k0-1)/2);


if isfield(ops,'scale_sigma')
    PD.scale_sigma = ops.scale_sigma;            
else
    PD.scale_sigma = 1/Normb;      
end

% this rescales c 
if isfield(ops,'scale_rho')
    PD.scale_rho = ops.scale_rho;            
else
    PD.scale_rho = 1;      
end

%% Anderson acceleration
% there are type 1 and type 2 options see SCS paper
% lookback determins how many past iterations are included
if isfield(ops,'AAlookback')
    PD.AAlookback = ops.AAlookback;            
else
    PD.AAlookback = 10;      
end

if isfield(ops,'AAtyp')
    PD.AAtyp = ops.AAtyp;            
else
    PD.AAtyp = 1;      
end

% this is a regularization parameter
if isfield(ops,'AAreg')
    PD.AAreg = ops.AAreg;            
else
    switch PD.AAtyp
        case 2
            PD.AAreg = 1e-12;
        case 1
            PD.AAreg = 1e-8;
    end
end

% another internal parameter
if isfield(ops,'AAsafeguard')
    PD.AAsafeguard = ops.AAsafeguard;            
else
    PD.AAsafeguard = 1;      
end



%% parallel pool settings
% this was just for experimenting with matlab parallel processing 
% there are multithread and multiprocess options. 
% multiprocess is of course useless for this kind of algorithm, as
% variables have to be sent to the different processes every time
% just use the default 'single' setting, which means no parallelization

if isfield(ops,'parallelize')
    PD.parallelize=ops.parallelize;
else
    PD.parallelize='single';
end

if isfield(ops,'numWorkers') % determines the number of workers in a PROCESS pool 
    PD.numWorkers = ops.numWorkers;
else
    PD.numWorkers = maxNumCompThreads;
end

if isfield(ops,'numThreads') % determines the number of threads in a THREADS pool *AND* the number of local theads
                             % on each worker in a PROCESS pool
    PD.numThreads = ops.numThreads;
else
    switch PD.parallelize
        case 'process'
            PD.numThreads = floor(maxNumCompThreads / PD.numWorkers);
        otherwise
            PD.numThreads = maxNumCompThreads;
    end
end

ps = parallel.Settings;
ps.Pool.AutoCreate = false;
prevPool=gcp
delete(prevPool)
PD.pool=[];

switch PD.parallelize
        case 'thread'
            if   PD.numThreads < maxNumCompThreads
                warning('Thread based pool will start with less threads than maxNumCompThreads!')
            end
            prevNumThreads = maxNumCompThreads(PD.numThreads); % reduce num threads
            PD.pool=parpool("threads");                        % open pool (number of threads = S.numThreads)
            maxNumCompThreads(prevNumThreads);                       % back to original max threads
            fprintf('Set max num comp threads back to %d \n',maxNumCompThreads)
        case 'process'
            PD.pool=parpool(PD.numWorkers);
            prevNumThreads = maxNumCompThreads;
            %set number of theads on all workers AND client
            cmd=['pctRunOnAll maxNumCompThreads(',num2str(PD.numThreads),')'];
            disp('setting maxNumCompThreads on all workers and cliens. Previous values were:')
            eval(cmd);
            disp('current values are:')
            pctRunOnAll maxNumCompThreads
            disp(['resetting client maxNumThreads to ',num2str(prevNumThreads)])
            maxNumCompThreads(prevNumThreads); %reset maxThreads on client to original value
        case 'single'
            disp('NO parpool started, to paralellize specify ops.parallelize="thread" or "process"')
        otherwise
            error('ops.parallelize parameter should be eiter process,thread or single')
end            
 



%% save progress warm start
% this allows to save the last iteration and then start another run from
% there 

if isfield(ops,'saveProgress')
    PD.saveProgress = ops.saveProgress;            
else
    PD.saveProgress = true;      
end

if  PD.saveProgress 
PD.filename=['OUT_SCSMF_Job',getenv('SLURM_JOB_ID'),'_',modelInfo.modelName,'_D',num2str(D),'_n',num2str(n),'_', datestr(now,30),'.mat'];
end

if isfield(ops,'saveStride')
    PD.saveStride = ops.saveStride;            
else
    PD.saveStride = 500;      
end

    
if isfield(ops,'warmStart')
    PD.warmStart = ops.warmStart;            
else
    % no field means don't warm start   
end

    


    
%% indicds for partial tracing 
% every iteration of SCS involves applying the linear constraints which
% involves partial traces of the matrix variables. This operation just sums
% up a bunch of entries of the respective variable. those indices are
% percomputed once


% for d=2   
% inds_L(:,1)=1:prod(dims)/2;
% inds_L(:,2)=(prod(dims)/2 + 1):prod(dims);
% inds_R(:,1)=1:2:prod(dims);
% inds_R(:,2)=inds_R(:,1)+1;

dims0 = d *ones(1,k0+1);
dims  = [d D D d];

PD.inds_L0 = reshape(1:prod(dims0),prod(dims0)/d,d);
PD.inds_R0 = reshape(1:prod(dims0),d,prod(dims0)/d)';
PD.inds_L = reshape(1:prod(dims),prod(dims)/d,d);
PD.inds_R = reshape(1:prod(dims),d,prod(dims)/d)';


%% save extended isometries 
% this saves computing the tensor products at each iteration

PD.V0xI=kron(V0,speye(d)); 
PD.IxV0=kron(speye(d),V0); 
PD.LxI=cell(1,n-1);
PD.IxR=cell(1,n-1);

for k=k0:n-1
    PD.LxI{k} = kron(L{k},speye(d));
    PD.IxR{k} = kron(speye(d),R{k});
end

%% save vars on workers
% irrelevant...
% in the case of parallel processing some things can be saved on each
% worker 
if strcmp(PD.parallelize,'process')
    % save LxI and IxR as constant vars on workers
    disp('saving L and R maps on workers')
    PD.LxIconst = parallel.pool.Constant(PD.LxI);
    PD.IxRconst = parallel.pool.Constant(PD.IxR);
    
    disp('saving partial trace indices on workers')
    PD.inds_Lconst = parallel.pool.Constant(PD.inds_L);
    PD.inds_Rconst = parallel.pool.Constant(PD.inds_R );
    
end

 
%% ground state energy of the model and the MPS upper bound
PD.Eexact=Eexact;
PD.Eupper=Eupper;


    