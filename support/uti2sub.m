function [i,j]=uti2sub(len,uti)
%UTI2SUB    Square matrix upper triangle linear index to subscripts

uti=uti(:);
k=cumsum([0 1:len-1]);
[i,j]=min((uti(:,ones(1,len))>k(ones(length(uti),1),:)).');
i=uti-k(j-1).'; j=j(:);

end
