function [S] = retrieveDataFromOUTfile(model,FolderPath)
if nargin < 2
    FolderPath = 'C:\Users\ilya\Dropbox\IQOQI\MarginalEnergy\TIspinChain\output\clusterRuns\SCSMF\';
end
dirContent=dir(FolderPath);
dirContent=dirContent(~[dirContent.isdir]);
dirContent=dirContent(endsWith({dirContent.name},'.out'));

if ~strcmp(model,'all')
    inds = find(contains({dirContent.name},['_',model,'_']));
else 
    inds =1:size(dirContent,1);
end

numFiles=length(inds);
S=struct('jobNo',[],'D',[],'n',[],'deltaE',[],'T_Affine',[],'maxIter',[]);

for j=1:length(inds)
    fileName = dirContent(inds(j)).name;
    str =fileread([FolderPath,fileName]);
    if contains(str,'Lower bound energy difference :') % means run finished
        C = regexp(str,'(?<=Lower bound energy difference : )\d.\d*','match'); 
        if ~isempty(C)
            S(j).deltaE= str2num(C{1}); 
            C = regexp(str,'(?(?<!CGmaxIter: )(?<=maxIter: )\d*)','match'); 
            S(j).maxIterTop = str2num(C{1});
            C=splitlines(str); 
            lineNum=find( contains(C,'Iteration ended'))-1;
            extractLine=split(C{lineNum},' ');
            S(j).maxIter = str2num(extractLine{1}); 
            
            if mod(S(j).maxIter,100), S(j).maxIter = S(j).maxIterTop; end
            C = regexp(str,'(?<=D: )\d*','match'); 
            S(j).D = str2num(C{1});
            
            C = regexp(str,'(?<=n: )\d*','match'); 
            S(j).n = str2num(C{1});
            
            S(j).jobNo = str2num(fileName(end-10:end-4));
            
            C = regexp(str,'(?<=: T_Aff = )(\d*[.]\d*|\d*)','match');
            S(j).T_Affine= str2num( C{1});
            
            C = regexp(fileName,'(?<=AAtyp)\d*','match');
            if ~isempty(C)
                S(j).AAtyp = str2num(C{1});

                C = regexp(fileName,'(?<=AAlb)\d*','match');
                S(j).AAlb = str2num(C{1});

                C = regexp(fileName,'(?<=AAsg)\d*','match');
                S(j).AAsg = str2num(C{1});
            end
            linesStr=splitlines(str);
            lineBelowLastIter = find(contains(linesStr,'Iteration ended'));
            lastIterLine=linesStr{lineBelowLastIter-1};
            C = regexp(lastIterLine,'[\s,0-9,.,\-,e]*\|','match');
            S(j).resid = str2num( [C{2}(1:end-1),C{3}(1:end-1),C{4}(1:end-1)]);
             
            C = regexp(str,'(?<=filename: )\''.*\''','match'); 
            S(j).matFile = C{1}(2:end-1);
        end
    end
end
S=S(~cellfun(@isempty,{S.deltaE}));
[~,I] = sort([S(:).D]);
S=S(I);