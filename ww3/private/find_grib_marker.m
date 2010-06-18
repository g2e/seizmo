%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scan for GRiB Marker (GRIB)                       %
% Added Jun 2002 to fix an issue with certain       %
% REANALYSIS GRiB files.                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iret=find_grib_marker(fid,ScreenDiag)

current_position=ftell(fid);

nope=1;
iret=0;

while nope
   gribmarker=fread(fid,4);
   if ~isempty(gribmarker)&length(gribmarker)==4
      if strcmp(char(gribmarker'),'GRIB')
	 % GRIB marker found

	 if ScreenDiag
     	    str=sprintf('   GRiB Record marker "GRIB" found at chars %d-%d',...
     			ftell(fid)-3,ftell(fid));
     	    disp(str);
	 end
	 iret=ftell(fid);
	 break
      else
	 fseek(fid,-3,0);
      end 
   else
      % eof reached
      iret=0;
      nope=0;
      break
   end
end 
