function [] = printTableLine(fieldsData,wNp)
% wNp := width.precision
str=[];
for i=1:numel(fieldsData)
    str = [str,'%-',num2str(wNp(i)),'g | '];
end
str=[str,'\n'];
fprintf(str,fieldsData)