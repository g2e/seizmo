function [mod]=readcrust10(bndsfile,rhofile,vpfile,vsfile)
%READCRUST10    Reads Crust1.0 files, putting them into a struct
%
% Usage: mod=readcrust10(...
%            'crust1.bnds','crust1.rho','crust1.vp','crust1.vs');
%
% Description:
%  mod.top is the top of each layer in km from sea level
%     - mod.top(:,:,1) gives the elevation/sealevel
%     - mod.top(:,:,2) gives the elevation/bathymetry
%     - mod.top(:,:,3) gives the elevation/bathymetry under ice
%     - mod.top(:,:,9) gives the moho depth
%     - sed thick   = 3:5
%     - xtl thick   = 6:8
%     - crust thick = 2:8
%  mod.rho is the density in kg/m3
%  mod.vp is p velocity in km/s
%  mod.vs is s velocity in km/s


% read in model putting it into a useful format
mod.top=permute(reshape(str2num(readtxt(bndsfile)),[360 180 9]),[2 1 3]);
mod.rho=permute(reshape(str2num(readtxt(rhofile)),[360 180 9]),[2 1 3]);
mod.vp=permute(reshape(str2num(readtxt(vpfile)),[360 180 9]),[2 1 3]);
mod.vs=permute(reshape(str2num(readtxt(vsfile)),[360 180 9]),[2 1 3]);

% FOR EASY ACCESS & STORAGE: make flat table & convert to integers
%mod.top=int16(reshape(mod.top,[360*180 9])*100);
%mod.rho=int16(reshape(mod.rho,[360*180 9])*100);
%mod.vp=int16(reshape(mod.vp,[360*180 9])*100);
%mod.vs=int16(reshape(mod.vs,[360*180 9])*100);

end
