%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test to see if fid points to a real GRiB File     %
% Added Jun 2002 to fix an issue with certain       %
% REANALYSIS GRiB files.                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iret=is_grib_file(fid,ScreenDiag)
nope=1;
iret=0;
if ScreenDiag
   disp('Searching for first GRIB marker...')
end

iret=find_grib_marker(fid,ScreenDiag);
if iret>0,iret=1;,end


