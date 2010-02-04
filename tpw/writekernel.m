function []=writekernel(filename,Kph,Kam,x,y)
% write out sensitivity kernel

% combine (note the switch)
m=[x(:) y(:) Kph(:) Kam(:)];

% write out header portion
n=size(x,1); d=y(2)-y(1);
dlmwrite(filename,[n x(1) d; n x(1) d],' ')

% write kernels
dlmwrite(filename,m,'-append','delimiter',' ');

end