%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLE 6, GDS (Section 2) Octet 6                 %
% BOB 28 Oct 2005 Added more DRT's                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dat_Rep_Type=table6(ival)
switch ival
   case 0,    Dat_Rep_Type = 'Equidis. Cyl. Lat/Lon';
   case 1,    Dat_Rep_Type = 'Mercator';
   case 2,    Dat_Rep_Type = 'Gnomic';
   case 3,    Dat_Rep_Type = 'Lambert Conf.';
   case 4,    Dat_Rep_Type = 'Gaussian Lat/Lon';
   case 5,    Dat_Rep_Type = 'Polar Stereogrphic';
   case 6,    Dat_Rep_Type = 'Universal Transverse Mercator (UTM) projection';
   case 7,    Dat_Rep_Type = 'Simple polyconic projection';
   case 8,    Dat_Rep_Type = 'Albers equal-area, secant or tangent, conic or bi-polar, projection';
   case 9,    Dat_Rep_Type = 'Millers cylindrical projection';
   case 10,   Dat_Rep_Type = 'Rotated latitude/longitude grid';
   case 13,   Dat_Rep_Type = 'Oblique Lambert Conf.';
   case 14,   Dat_Rep_Type = 'Rotated Gaussian latitude/longitude grid';
   case 50,   Dat_Rep_Type = 'Spher. Harm. Coeff';
   case 90,   Dat_Rep_Type = 'Space View';
   case 201,  Dat_Rep_Type = 'Arakawa Semi-Staggered E-Grid';
   case 202,  Dat_Rep_Type = 'Arakawa Filled E-Grid';
   otherwise, Dat_Rep_Type = 'Reserved';
end
