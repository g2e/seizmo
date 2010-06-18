%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLE 0, PDS Octet 5                             %
% National/International Originating Centers       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function orig_center=table0(ival)
switch ival
   case  1,   orig_center='Melbourne (WMC)';
   case  2,   orig_center='Melbourne (WMC)';
   case  4,   orig_center='Moscow (WMC)';
   case  5,   orig_center='Moscow (WMC)';
   case  7,   orig_center='NCEP (WMC)';
   case  8,   orig_center='NWSTG (WMC)';
   case  9,   orig_center='Other (WMC)';
   case 34,   orig_center='Japan. Met. Soc. (RSMC)';
   case 52,   orig_center='Nat. Hurr. Center, Miami';
   case 53,   orig_center='Can. Meteor. Serv.-Montreal (RSMC)';
   case 57,   orig_center='US Air Force-Global Weather Center';
   case 58,   orig_center='US Navy-FNOC';
   case 59,   orig_center='NOAA Forecast Systems Lab. (Boulder)';
   case 60,   orig_center='NCAR, Boulder';
   case 74,   orig_center='U.K. Met Office - Bracknell';
   case 85,   orig_center='French Weather Service - Toulouse';
   case 97,   orig_center='European Space Agency (ESA)';
   case 98,   orig_center='ECMWF - Reading';
   case 99,   orig_center='DeBilt, Netherlands';
   otherwise, orig_center='N/A';
end
