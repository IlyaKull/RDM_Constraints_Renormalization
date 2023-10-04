Ds=6:12;
ns=[50 100 200 400];
ts=cell(length(Ds),length(ns));
for j=1:length(Ds)
    for i=1:length(ns)
        ts{j,i}=testApplyConstraints(1,Ds(j),ns(i));
    end
end

[time,winner]=cellfun(@min,ts,'UniformOutput',false)
ratio=@(x) x(2)/x(1);
rat=cellfun(ratio,ts,'UniformOutput',false)