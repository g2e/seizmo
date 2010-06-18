%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLE 4, PDS (Section 1) Octet 18                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function time_unit=table4(ival)
switch ival
   case  0,time_unit='minute';
   case  1, time_unit='hour';
   otherwise,time_unit='none';
end
