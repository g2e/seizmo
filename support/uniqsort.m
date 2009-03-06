function [x]=uniqsort(x)
%UNIQSORT    A fast unique sort for numeric vectors

x=sort(x(:));
x([diff(x(:)); true]==0)=[];

end