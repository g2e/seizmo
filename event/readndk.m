function [events]=readndk(file)
%READNDK    Reads a GlobalCMT Project NDK-format text file into a struct
%
%    Usage:    events=readndk(file)
%
%    Description: EVENTS=READNDK(FILE) reads in an NDK-formatted text file
%     from the Global CMT project (www.globalcmt.org).  All event info from
%     the file is imported into struct array EVENTS (see the Notes section
%     below for more details).  If FILE is not given or set to '' then a
%     graphical file selection menu is presented.
%
%    Notes:
%     - Details of the NDK format may be found using the Global CMT
%       project's website (www.globalcmt.org) or (hopefully) here:
%       http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/
%       The file 'allorder.ndk_explained' will provide the details.
%     - Fields of NDK struct:
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
%           strik2
%           dip2
%           rake2
%
%    Examples:
%     Import info from a quick CMT into some records:
%      ndk=readndk('quick.ndk');
%      data=setevent(data,ndk(33));
%
%    See also: READSODEVENTCSV, READTXT, SETEVENT

%     Version History:
%        Mar.  6, 2009 - initial version (in SEIZMO)
%        Apr. 23, 2009 - octave compatibility fix
%        Jan. 26, 2010 - allow no input (select file graphically), added
%                        history and documentation, clean up code a bit
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2010 at 22:40 GMT

% todo:

% check nargin
msg=nargchk(0,1,nargin);
if(~isempty(msg)); error(msg); end;

% graphical selection
if(nargin<1 || isempty(file))
    [file,path]=uigetfile(...
        {'*.ndk;*.NDK' 'NDK Files (*.ndk,*.NDK)';
        '*.*' 'All Files (*.*)'},...
        'Select NDK File');
    if(isequal(0,file))
        error('seizmo:readndk:noFileSelected','No input file selected!');
    end
    file=strcat(path,filesep,file);
end

% check file
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

% open file for reading
fid=fopen(file);

% check if file is openable
if(fid<0)
    error('seizmo:readndk:cannotOpenFile',...
        'NDK File: %s\nNot Openable!',file);
end

% read in file and close
a=fread(fid,'*char');
fclose(fid);

% delete linefeed and carriage return characters
a(a==10)=[];
a(a==13)=[];

try
    % reshape into nx80 array
    a=reshape(a,80,[])';
catch
    error('seizmo:readndk:malformedNDK',...
        ['NDK file: %s\n'...
        'Looks like some lines are not 80 characters long!'],file);
end

% check that there are 5 lines per event
if(mod(size(a,1),5))
    error('seizmo:readndk:malformedNDK',...
        ['NDK file: %s\n'...
        'Looks like some events do not have 5 lines!'],file);
end

% allocate struct
n=size(a,1)/5;
events(1:n,1)=struct();

% begin slow section (text to datatype to struct field)
% - replacement with a faster version is welcome!
% LINE 1
temp=cellstr(a(1:5:end,1:4));
[events.catalog]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,6:9)));
[events.year]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,11:12)));
[events.month]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,14:15)));
[events.day]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,17:18)));
[events.hour]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,20:21)));
[events.minute]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,23:26)));
[events.seconds]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,28:33)));
[events.latitude]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,35:41)));
[events.longitude]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,43:47)));
[events.depth]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,49:51)));
[events.mb]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,53:55)));
[events.ms]=deal(temp{:});
temp=cellstr(a(1:5:end,57:80));
[events.location]=deal(temp{:});
% LINE 2
temp=cellstr(a(2:5:end,1:16));
[events.name]=deal(temp{:});
temp=cellstr(a(2:5:end,18));
[events.data1type]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,20:22)));
[events.data1nstn]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,23:27)));
[events.data1ncmp]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,28:31)));
[events.data1sper]=deal(temp{:});
temp=cellstr(a(2:5:end,33));
[events.data2type]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,35:37)));
[events.data2nstn]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,38:42)));
[events.data2ncmp]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,43:46)));
[events.data2sper]=deal(temp{:});
temp=cellstr(a(2:5:end,48));
[events.data3type]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,50:52)));
[events.data3nstn]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,53:57)));
[events.data3ncmp]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,58:61)));
[events.data3sper]=deal(temp{:});
temp=cellstr(a(2:5:end,63:68));
[events.srcinvtype]=deal(temp{:});
temp=cellstr(a(2:5:end,70:74));
[events.srcfunctype]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,76:80)));
[events.srcfuncdur]=deal(temp{:});
% LINE 3
temp=cellstr(a(3:5:end,1:8));
[events.centroidstring]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,13:18)));
[events.centroidtime]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,20:22)));
[events.centroidtimeerr]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,24:29)));
[events.centroidlat]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,31:34)));
[events.centroidlaterr]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,36:42)));
[events.centroidlon]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,44:47)));
[events.centroidlonerr]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,49:53)));
[events.centroiddep]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,55:58)));
[events.centroiddeperr]=deal(temp{:});
temp=cellstr(a(3:5:end,60:63));
[events.depthtype]=deal(temp{:});
temp=cellstr(a(3:5:end,65:80));
[events.timestamp]=deal(temp{:});
% LINE 4
temp=num2cell(str2num(a(4:5:end,1:2)));
[events.exponent]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,4:9)));
[events.mrr]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,11:15)));
[events.mrrerr]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,17:22)));
[events.mtt]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,24:28)));
[events.mtterr]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,30:35)));
[events.mpp]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,37:41)));
[events.mpperr]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,43:48)));
[events.mrt]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,50:54)));
[events.mrterr]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,56:61)));
[events.mrp]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,63:67)));
[events.mrperr]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,69:74)));
[events.mtp]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,76:80)));
[events.mtperr]=deal(temp{:});
% LINE 5
temp=cellstr(a(5:5:end,1:3));
[events.version]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,5:11)));
[events.eigval1]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,13:14)));
[events.plunge1]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,16:18)));
[events.azimuth1]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,20:26)));
[events.eigval2]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,28:29)));
[events.plunge2]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,31:33)));
[events.azimuth2]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,35:41)));
[events.eigval3]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,43:44)));
[events.plunge3]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,46:48)));
[events.azimuth3]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,50:56)));
[events.scalarmoment]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,58:60)));
[events.strike1]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,62:63)));
[events.dip1]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,65:68)));
[events.rake1]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,70:72)));
[events.strik2]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,74:75)));
[events.dip2]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,77:80)));
[events.rake2]=deal(temp{:});

end
