function [mod]=readcrust10(bndsfile,rhofile,vpfile,vsfile)
%READCRUST10    Reads Crust1.0 files, putting them into a struct
%
%    Usage:    mod=readcrust10(bnds,rho,vp,vs)
%
%    Description:
%     MOD=READCRUST10(BNDS,RHO,VP,VS) reads in the Crust1.0 model from the
%     files given by BNDS, RHO, VP & VS.  Note that the output is a struct
%     with the model values stored in integer format for memory savings.
%     This format needs to be converted to double precision and divided by
%     100 to get the actual values.  Please use GETCRUST to access the
%     values as this is not really meant to be used for routine access.
%
%    Notes:
%
%    Examples:
%     % Read in the model from an unzipped Crust1.0 archive:
%     mod=readcrust10('crust1.bnds','crust1.rho','crust1.vp','crust1.vs');
%
%    See also: READCRUST2, GETCRUST

%     Version History:
%        Jan. 23, 2014 - initial version, added full docs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 23, 2014 at 02:45 GMT

% todo:

% read in model putting it into a useful format
mod.top=permute(reshape(str2num(readtxt(bndsfile)),[360 180 9]),[2 1 3]);
mod.rho=permute(reshape(str2num(readtxt(rhofile)),[360 180 9]),[2 1 3]);
mod.vp=permute(reshape(str2num(readtxt(vpfile)),[360 180 9]),[2 1 3]);
mod.vs=permute(reshape(str2num(readtxt(vsfile)),[360 180 9]),[2 1 3]);

% FOR EASY ACCESS & STORAGE: make flat table again & convert to integers
mod.top=int16(reshape(mod.top,[360*180 9])*100);
mod.rho=int16(reshape(mod.rho,[360*180 9])*100);
mod.vp=int16(reshape(mod.vp,[360*180 9])*100);
mod.vs=int16(reshape(mod.vs,[360*180 9])*100);

end
