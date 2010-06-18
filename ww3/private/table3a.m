%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLE 3a, PDS Octet 10  (Special Levels)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ret1,ret2]=table3a(ival)
ret2='';
switch ival
   case 1,   ret1='Surface';                            ret2='SFC';
   case 2,   ret1='Cloud Base Level';                   ret2='CBL';
   case 3,   ret1='Cloud Top Level';                    ret2='CTL';
   case 4,   ret1='Level of 0 deg (C) isotherm';        ret2='0DEG';
   case 5,   ret1='Level of Adiabatic Condensation';    ret2='ADCL';
   case 6,   ret1='Maximum wind level';                 ret2='MWSL';
   case 7,   ret1='Tropopause';                         ret2='TRO';
   case 8,   ret1='Nominal Top of Atmos';               ret2='NTAT';
   case 9,   ret1='Sea Bottom';                         ret2='SEAB';
   case 20,  ret1='Isothermal level';                   ret2='TMPL';
   case 100, ret1='Isobaric surface Pressure in hPa (2 octets)'; 
   case 101, ret1={'Layer between two isobaric surfaces','Pressure of top in kPa','Pressure of bottom in kPa'};
   case 102, ret1='Mean sea level 0';
   case 103, ret1={'Specified altitude above mean sea level','Altitude in meters (2 octets)'};
   case 104, ret1={'Layer between two specified altitudes above mean sea level','Altitude of top in hm','Altitude of bottom in hm'};
   case 105, ret1={'Specified height above ground','Height in meters (2 octets)'};
   case 106, ret1={'Layer between two specified height levels above ground','Height of top in hm','Height of bottom in hm'};
   case 107, ret1={'Sigma level','Sigma value in 1/10000 (2 octets)'};
   case 108, ret1={'Layer between two sigma levels','Sigma value of top in 1/100','Sigma value of bottom in 1/100'};
   case 109, ret1={'Hybrid level','Level number (2 octets)'};
   case 110, ret1={'Layer between two hybrid levels','Level number of top','Level number of bottom'};
   case 111, ret1={'Departmenth below land surface','Departmenth in centimeters (2 octets)'};
   % NCEP Special levels and Layers'
   case 204, ret1='Highest tropospheric freezing level';ret2='HTFL';
   case 209, ret1='Boundary layer cloud bottom level '; ret2='BCBL';
   case 210, ret1='Boundary layer cloud top level';     ret2='BCTL';
   case 211, ret1='Boundary layer cloud layer';         ret2='BCY';
   case 212, ret1='Low cloud bottom level';             ret2='LCBL';
   case 213, ret1='Low cloud top level';                ret2='LCTL';
   case 214, ret1='Low cloud layer';                    ret2='LCY';
   case 222, ret1='Middle cloud bottom level';          ret2='MCBL';
   case 223, ret1='Middle cloud top level';             ret2='MCTL';
   case 224, ret1='Middle cloud layer';                 ret2='MCY';
   case 232, ret1='High cloud bottom level';            ret2='HCBL';
   case 233, ret1='High cloud top level';               ret2='HCTL';
   case 234, ret1='High cloud layer';                   ret2='HCY';
   case 242, ret1='Convective cloud bottom level';      ret2='CCBL';
   case 243, ret1='Convective cloud top level';         ret2='CCTL';
   case 244, ret1='Convective cloud layer';             ret2='CCY';
   otherwise, ret1='Reserved';ret2='';   
end
