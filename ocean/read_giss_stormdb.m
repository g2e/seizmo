function [storms]=read_giss_stormdb(file)
%READ_GISS_STORMDB    Reads GISS Atlas of Extratropical Storm Tracks
%
%    Usage:    storms=read_giss_stormdb(file)
%
%    Description:
%     STORMS=READ_GISS_STORMDB(FILE) reads in the storm track data from
%     NASA's Goddard Institute for Space Studies Atlas of Extratropical
%     Storm Tracks.  The storm tracks are contained in a binary file with
%     unformatted sequential records.  The output STORMS is a struct with
%     the following fields:
%       STORMS.time     -- times of track points (in datenum format)
%             .lat      -- latitudes of track points (in degrees)
%             .lon      -- longitudes of track points (in degrees)
%             .pressure -- pressure at track points (in mbars)
%     Each index of STORMS contains the track for a separate storm so
%     STORMS(3) contains the track data for the third storm in the
%     database.
%
%    Notes:
%     - File is available at:
%        http://data.giss.nasa.gov/stormtracks/stormtracks.bin
%
%    Examples:
%     % You may omit the file input to graphically select the file:
%     storms=read_giss_stormdb;
%
%    See also: READ_HURDAT, MAPSTORMS

%     Version History:
%        Feb. 12, 2013 - initial version
%        Feb. 16, 2013 - storms now returned as scalar struct
%        Jan. 26, 2014 - abs path exist fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2014 at 13:30 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% filterspec
filterspec={'*.bin;*.BIN' 'BIN Files (*.bin,*.BIN)';
    '*.*' 'All Files (*.*)'};

% graphical selection
if(nargin<1 || isempty(file))
    [file,path]=uigetfile(filterspec,'Select File');
    if(isequal(0,file))
        error('seizmo:read_giss_stormdb:noFileSelected',...
            'No input file selected!');
    end
    file=strcat(path,filesep,file);
else
    % check file
    if(~isstring(file))
        error('seizmo:read_giss_stormdb:fileNotString',...
            'FILE must be a string!');
    end
    if(~isabspath(file)); file=[pwd fs file]; end
    if(~exist(file,'file'))
        error('seizmo:read_giss_stormdb:fileDoesNotExist',...
            'File: %s\nDoes Not Exist!',file);
    elseif(exist(file,'dir'))
        error('seizmo:read_giss_stormdb:dirConflict',...
            'File: %s\nIs A Directory!',file);
    end
end

% open file (big endian)
fid=fopen(file,'rb','b');

% file size
fseek(fid,0,'eof');
bytes=ftell(fid);
fseek(fid,0,'bof');

% preallocate struct (guessing 50000 storms)
guess=50000;
s(guess,1)=struct('time',[],'lat',[],'lon',[],'pressure',[]);

% loop over times entries until done
j=0;
%print_time_left(0,bytes);
while(ftell(fid)<bytes)
    fseek(fid,4,'cof'); % skip record #
    x=fread(fid,2,'*int32'); % time (YYYYMMDDHH), # of storms
    fseek(fid,4,'cof'); % skip record #
    %print_time_left(ftell(fid),bytes);
    if(~x(2)); continue; end % any storms?
    for i=1:x(2) % loop over storms
        fseek(fid,4,'cof'); % skip record #
        nt=fread(fid,1,'int32'); % # of times for storm
        y=fread(fid,3*nt,'float32'); % lat (deg),lon (deg),pressure (mbar)
        fseek(fid,4,'cof'); % skip record #
        %print_time_left(ftell(fid),bytes);
        
        % get times
        j=j+1;
        s(j).time=(0:nt-1)'/2+x2time(x(1));
        
        % lat, lon, pressure
        s(j).lat=y(1:3:end);
        s(j).lon=y(2:3:end);
        s(j).pressure=y(3:3:end);
    end
end

% close file
fclose(fid);

% trim excess if any
if(j<guess); s(j+1:end)=[]; end

% compress into a scalar struct
storms.time={s.time}';
storms.lat={s.lat}';
storms.lon={s.lon}';
storms.pressure={s.pressure}';

end


function [t]=x2time(x)
% convert YYYYMMDDHH to datenum number
t=zeros(1,6,'int32');
ymd=x/100;
t(4)=x-ymd*100;
ym=ymd/100;
t(3)=ymd-ym*100;
t(1)=ym/100;
t(2)=ym-t(1)*100;
t=datenum(double(t));
end
