%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% int3                                                      %
%                                                           %
% BOB 6 Sept 2002:  MATLAB6.5 (R13) has an improved         %
% HEX2DEC function that actually forces the input to        %
% be a string. Hence, the call hex2dec(80) causes an        %
% error in R13. The lack of an error in previous            %
% versions was fortuitous since the (unnecessary)           %
% logic in the lines                                        %
%                                                           %
% sgn=1;                                                    %
% tmp=dec2bin(a,8);if strcmp(tmp(1),'1'),sgn=-1;,end        %
%                                                           %
% essentially negated the effect of the hex2dec(80) call.   %
%                                                           %
% Therefore, the line                                       %
%  retval=sgn*((1-(bitshift(bitand(a,hex2dec(80)),-6)))*... %
%  (bitshift(bitand(a,127),16)+(bitshift(b,8)+c)));         %
%                                                           %
% has been replaced with the line below, since              %
% hex2dec('80')=128                                         %
%                                                           %
% The same problem exists in the routine int2 below,        %
% and is similarly fixed.                                   %
%                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval=int3(a,b,c)     
retval=(1-(bitshift(bitand(a,128),-6)))*...
          (bitshift(bitand(a,127),16)+...
          (bitshift(b,8)+c));    
