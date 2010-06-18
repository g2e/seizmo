%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLE 5, PDS (Section 1) Octet 21                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tri_string=table5(ival)
switch ival
   case  0,tri_string='Valid at P1';
   case  1,tri_string='Anal Prod for P1';
   case  2,tri_string='P1<Valid<P2';
   case  3,tri_string='Average';
   case  4,tri_string='Accumulation';
   case  5,tri_string='Difference';
   case  51,tri_string='Clim. Mean Value';
   otherwise,tri_string='N/A';
end
