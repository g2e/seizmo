function []=writekernel(filename,x,y,sp,sa)
% write out sensitivity kernel

% combine (note the switch)
m=[y(:) x(:) sp(:) sa(:)];

% write out header portion
n=size(x,1);
dx=x(2)-x(1);
dlmwrite(filename,[n x(1) dx; n x(1) dx],' ')

% write kernels
dlmwrite(filename,m,'-append','delimiter',' ');

end