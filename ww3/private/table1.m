%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLE 1, PDS Octet 8                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ret=table1(ival)
bval=dec2bin(ival);
if strcmp(bval(1),'0'),ret=0;,end
if strcmp(bval(1),'1'),ret=1;,end
if strcmp(bval(2),'0'),ret=[ret 0];,end
if strcmp(bval(2),'1'),ret=[ret 1];,end
