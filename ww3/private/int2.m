%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% int2                                              %
%                                                   %
% See notes in int3.m                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval=int2(a,b)                                                      
retval=((1-(bitshift(bitand(a,128),-6)))*...
           (bitshift(bitand(a,127),8)+b)); 
