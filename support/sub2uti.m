function [uti]=sub2uti(i,j)
%SUB2UTI    Square matrix upper triangle linear index from subscripts

j=j-1;
k=cumsum([0 1:max(j)]);
uti=k(j).'+i(:);

end
