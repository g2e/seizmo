function [i,j]=lti2sub(len,lti)
%LTI2SUB    Square matrix lower triangle linear index to subscripts

lti=lti(:);
k=cumsum([0 len-1:-1:1]);
[i,j]=min((lti(:,ones(1,len))>k(ones(length(lti),1),:)).');
j=j(:)-1; k=cumsum(1:len-1); i=lti+k(j).'-len*(j-1);

end
