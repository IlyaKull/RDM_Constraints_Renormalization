function []= makeTableHeader(fields,width)
str=[];
for i=1:numel(fields)
    str=[str,pad(fields{i},width),' | '];
end
fprintf([str,'\n'])
