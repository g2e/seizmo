%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% EXTRACT GRiB RECORD %%%%%%%%%%%%%%%
%
% 08 Feb 2007 : added check for thinned/reduced GG's
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pds_struct,gds_struct,bms_struct,bds_struct,dataarray]=...
         extract_grib(fid,nrec,fpos,headerskip,dataskip)

% read section 1, the PDS (Product Definition Section)
oct1to3=fread(fid,3);
lenpds=bitshift3(oct1to3(1),oct1to3(2),oct1to3(3));
fseek(fid,-3,0);
pds_struct=get_pds(fid,lenpds);

% read section 2, the GDS (Grid Definition Section)
if pds_struct.HASGDS
   oct1to3=fread(fid,3);
   lengds=bitshift3(oct1to3(1),oct1to3(2),oct1to3(3));
   fseek(fid,-3,0);
   gds_struct=get_gds(fid,lengds);
else
   gds_struct=[];
end

% read section 3, the BMS (Bit Map Section)
if pds_struct.HASBMS
   oct1to3=fread(fid,3);
   lenbms=bitshift3(oct1to3(1),oct1to3(2),oct1to3(3));
   fseek(fid,-3,0);
   bms_struct=get_bms(fid,lenbms);
else
   bms_struct.bitmap=uint8([]);
end

% read section 4, the MANDATORY BDS (Binary Data Section)
oct1to3=fread(fid,3);
lenbds=bitshift3(oct1to3(1),oct1to3(2),oct1to3(3));
fseek(fid,-3,0);
if ~dataskip
   bds_struct=get_bds(fid,lenbds);
else
   bds_struct=[];
   fseek(fid,lenbds,0);
end

% Its time to decode the binary data in bds_struct.bindata
if ~dataskip
   decimal_sf=10^(-pds_struct.DecSF);
   bin_sf=2^(bds_struct.bsfE);
   if isfield(gds_struct,'Ni')
      % 08 Feb 2007 : added check for thinned/reduced GG's
      if isfield(gds_struct,'NxNy')
         nxny=gds_struct.NxNy;
      else
         nxny=gds_struct.Ni*gds_struct.Nj; % Gaussian, etc, Grids
      end
   else
      nxny=gds_struct.Nx*gds_struct.Ny; % Lambert Grids, etc.
   end
   fac1=decimal_sf*bds_struct.RefVal;
   fac2=decimal_sf*bin_sf;
   dataarray=BDS_unpack_mex5(bds_struct.bindata,...
     			     bms_struct.bitmap,...
     			     bds_struct.nbits,...
     			     nxny,fac1,fac2);
else
   dataarray=[];
end
