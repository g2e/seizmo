function [lti]=sub2lti(len,i,j)
%SUB2LTI    Square matrix lower triangle linear index from subscripts

k=cumsum(1:len-1);
lti=i(:)-k(j).'+len*(j(:)-1);

end
