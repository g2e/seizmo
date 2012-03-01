function [events]=readndk(file,flag)
%READNDK    Reads a GlobalCMT Project NDK-format text file into a struct
%
%    Usage:    events=readndk(file)
%              events=readndk(string,true)
%
%    Description:
%     EVENTS=READNDK(FILE) reads in an NDK-formatted text file from the
%     Global CMT project (www.globalcmt.org).  All event info from the file
%     is imported into the struct EVENTS (see the Notes section below for
%     more details).  If FILE is not given or set to '' then a graphical
%     file selection menu is presented.
%
%     EVENTS=READNDK(STRING,TRUE) will take the first input to be a string
%     of an NDK file (rather than the filename) if the second input is set
%     to logical TRUE.  The string should be the same as if the NDK file
%     was read with READTXT (a single row char vector with linefeeds
%     included).
%
%    Notes:
%     - Details of the NDK format may be found using the Global CMT
%       project's website (www.globalcmt.org) or (hopefully) here:
%       http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/
%       The file 'allorder.ndk_explained' will provide the details.
%     - EVENTS is a scalar struct with many fields.  All fields are column
%       vectors (either of double or cellstr type) with as many rows as
%       events in the NDK file.  The fields of the NDK struct:
%           catalog
%           year
%           month
%           day
%           hour
%           minute
%           seconds
%           latitude
%           longitude
%           depth
%           mb
%           ms
%           location
%           name
%           data1type
%           data1nstn
%           data1ncmp
%           data1sper
%           data2type
%           data2nstn
%           data2ncmp
%           data2sper
%           data3type
%           data3nstn
%           data3ncmp
%           data3sper
%           srcinvtype
%           srcfunctype
%           srcfuncdur
%           centroidstring
%           centroidtime
%           centroidtimeerr
%           centroidlat
%           centroidlaterr
%           centroidlon
%           centroidlonerr
%           centroiddep
%           centroiddeperr
%           depthtype
%           timestamp
%           exponent
%           mrr
%           mrrerr
%           mtt
%           mtterr
%           mpp
%           mpperr
%           mrt
%           mrterr
%           mrp
%           mrperr
%           mtp
%           mtperr
%           version
%           eigval1
%           plunge1
%           azimuth1
%           eigval2
%           plunge2
%           azimuth2
%           eigval3
%           plunge3
%           azimuth3
%           scalarmoment
%           strike1
%           dip1
%           rake1
%           strike2
%           dip2
%           rake2
%
%    Examples:
%     % Read in the included example NDK file:
%     ndk=readndk('mar10.ndk');
%
%     % Extract a specific CMT (replace CMTIDX):
%     cmt=ssidx(ndk,CMTIDX);
%
%     % Import info from the CMT into some records:
%     data=setevent(data,cmt);
%
%    See also: READSODEVENTCSV, READTXT, SETEVENT, SSIDX, PARSE_ISC_ORIGIN,
%              FINDCMT, FINDCMTS, GLOBALCMT_UPDATE

%     Version History:
%        Mar.  6, 2009 - initial version (in SEIZMO)
%        Apr. 23, 2009 - octave compatibility fix
%        Jan. 26, 2010 - allow no input (select file graphically), added
%                        history and documentation, clean up code a bit
%        July 30, 2010 - change strik2 to strike2, now outputs a scalar
%                        struct, nargchk fix, use readtxt/getwords
%        Aug.  2, 2010 - updated catalog examples
%        Aug.  3, 2010 - allow string input
%        Mar.  1, 2012 - remove outdated examples
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  1, 2012 at 15:40 GMT

% todo:
% - Currently no numeric fields have missing values that would need to be
%   filled with nans.  If this does occur then use READNDK_SLOW method for
%   the field(s) with the missing value(s) but not for the rest as it is
%   quite slow.

% check nargin
error(nargchk(0,2,nargin));

% default flag
if(nargin<2 || isempty(flag)); flag=false; end

% skip if string input
if(~flag)
    % graphical selection
    if(nargin<1 || isempty(file))
        [file,path]=uigetfile(...
            {'*.ndk;*.NDK' 'NDK Files (*.ndk,*.NDK)';
            '*.*' 'All Files (*.*)'},...
            'Select NDK File');
        if(isequal(0,file))
            error('seizmo:readndk:noFileSelected',...
                'No input file selected!');
        end
        file=strcat(path,filesep,file);
    else % check file
        if(~ischar(file))
            error('seizmo:readndk:fileNotString',...
                'FILE must be a string!');
        end
        if(~exist(file,'file'))
            error('seizmo:readndk:fileDoesNotExist',...
                'CSV File: %s\nDoes Not Exist!',file);
        elseif(exist(file,'dir'))
            error('seizmo:readndk:dirConflict',...
                'CSV File: %s\nIs A Directory!',file);
        end
    end
    
    % read in ndk file
    txt=readtxt(file);
else
    % just copy file to txt
    if(nargin<1 || isempty(file))
        error('seizmo:readndk:emptyStr',...
            'STRING must be non-empty!');
    else
        txt=file;
    end
end

% delete carriage return characters
txt(txt==13)=[];

% split lines
lines=getwords(txt,sprintf('\n'));

% check that there are 80 char per line
if(any(cellfun('prodofsize',lines)~=80))
    error('seizmo:readndk:malformedNDK',...
        ['NDK file: %s\n'...
        'Looks like some lines are not 80 characters long!'],file);
end

% check that there are 5 lines per event
if(mod(numel(lines),5))
    error('seizmo:readndk:malformedNDK',...
        ['NDK file: %s\n'...
        'Looks like some events do not have 5 lines!'],file);
end

% now go back to a char array so we can extract certain column sections
txt=char(lines);

% begin slow section (text to datatype to struct field)
% - replacement with a faster version is welcome!
% LINE 1
events.catalog=cellstr(txt(1:5:end,1:4));
events.year=str2num(txt(1:5:end,6:9));
events.month=str2num(txt(1:5:end,11:12));
events.day=str2num(txt(1:5:end,14:15));
events.hour=str2num(txt(1:5:end,17:18));
events.minute=str2num(txt(1:5:end,20:21));
events.seconds=str2num(txt(1:5:end,23:26));
events.latitude=str2num(txt(1:5:end,28:33));
events.longitude=str2num(txt(1:5:end,35:41));
events.depth=str2num(txt(1:5:end,43:47));
events.mb=str2num(txt(1:5:end,49:51));
events.ms=str2num(txt(1:5:end,53:55));
events.location=cellstr(txt(1:5:end,57:80));
% LINE 2
events.name=cellstr(txt(2:5:end,1:16));
events.data1type=cellstr(txt(2:5:end,18));
events.data1nstn=str2num(txt(2:5:end,20:22));
events.data1ncmp=str2num(txt(2:5:end,23:27));
events.data1sper=str2num(txt(2:5:end,28:31));
events.data2type=cellstr(txt(2:5:end,33));
events.data2nstn=str2num(txt(2:5:end,35:37));
events.data2ncmp=str2num(txt(2:5:end,38:42));
events.data2sper=str2num(txt(2:5:end,43:46));
events.data3type=cellstr(txt(2:5:end,48));
events.data3nstn=str2num(txt(2:5:end,50:52));
events.data3ncmp=str2num(txt(2:5:end,53:57));
events.data3sper=str2num(txt(2:5:end,58:61));
events.srcinvtype=cellstr(txt(2:5:end,63:68));
events.srcfunctype=cellstr(txt(2:5:end,70:74));
events.srcfuncdur=str2num(txt(2:5:end,76:80));
% LINE 3
events.centroidstring=cellstr(txt(3:5:end,1:8));
events.centroidtime=str2num(txt(3:5:end,13:18));
events.centroidtimeerr=str2num(txt(3:5:end,20:22));
events.centroidlat=str2num(txt(3:5:end,24:29));
events.centroidlaterr=str2num(txt(3:5:end,31:34));
events.centroidlon=str2num(txt(3:5:end,36:42));
events.centroidlonerr=str2num(txt(3:5:end,44:47));
events.centroiddep=str2num(txt(3:5:end,49:53));
events.centroiddeperr=str2num(txt(3:5:end,55:58));
events.depthtype=cellstr(txt(3:5:end,60:63));
events.timestamp=cellstr(txt(3:5:end,65:80));
% LINE 4
events.exponent=str2num(txt(4:5:end,1:2));
events.mrr=str2num(txt(4:5:end,4:9));
events.mrrerr=str2num(txt(4:5:end,11:15));
events.mtt=str2num(txt(4:5:end,17:22));
events.mtterr=str2num(txt(4:5:end,24:28));
events.mpp=str2num(txt(4:5:end,30:35));
events.mpperr=str2num(txt(4:5:end,37:41));
events.mrt=str2num(txt(4:5:end,43:48));
events.mrterr=str2num(txt(4:5:end,50:54));
events.mrp=str2num(txt(4:5:end,56:61));
events.mrperr=str2num(txt(4:5:end,63:67));
events.mtp=str2num(txt(4:5:end,69:74));
events.mtperr=str2num(txt(4:5:end,76:80));
% LINE 5
events.version=cellstr(txt(5:5:end,1:3));
events.eigval1=str2num(txt(5:5:end,5:11));
events.plunge1=str2num(txt(5:5:end,13:14));
events.azimuth1=str2num(txt(5:5:end,16:18));
events.eigval2=str2num(txt(5:5:end,20:26));
events.plunge2=str2num(txt(5:5:end,28:29));
events.azimuth2=str2num(txt(5:5:end,31:33));
events.eigval3=str2num(txt(5:5:end,35:41));
events.plunge3=str2num(txt(5:5:end,43:44));
events.azimuth3=str2num(txt(5:5:end,46:48));
events.scalarmoment=str2num(txt(5:5:end,50:56));
events.strike1=str2num(txt(5:5:end,58:60));
events.dip1=str2num(txt(5:5:end,62:63));
events.rake1=str2num(txt(5:5:end,65:68));
events.strike2=str2num(txt(5:5:end,70:72));
events.dip2=str2num(txt(5:5:end,74:75));
events.rake2=str2num(txt(5:5:end,77:80));

end
