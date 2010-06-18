%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ibm2flt2                                          %
%                                                   %
% 30 Aug, 2005.                                     %
% this is a matlab version of the mex code          %
% ibm2fltmex5.c.  One less mex file to worry about. %
% code essentially from  wesley ebisuzaki, NCEP,    %
% wgrib, http://wesley.wwb.noaa.gov                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval=ibm2flt2(ibm)

mant = bitshift(ibm(2),16) + bitshift(ibm(3),8) + ibm(4);
if (mant == 0), retval=0.0;return;end

positive = (bitand(ibm(1), 128) == 0);

power =  bitand(ibm(1),127) - 64;

% c-version abspower = power > 0 ? power : -power;
abspower = abs(power);

expp=16.0;
retval=1.0;
while (abspower)
   if bitand(abspower,1)
      retval=retval*expp;
%      disp(sprintf('retval=%f',retval))
   end
   expp=expp*expp;
%   disp(sprintf('abspower=%d exp=%f retval=%f',abspower,expp,retval))
   abspower=bitshift(abspower,-1);
end

if (power < 0), retval = 1.0 / retval;,end
retval = retval * mant / 16777216.0;

if (positive == 0), retval = -retval;,end
%disp(sprintf('retval=%f\n',retval))

return
