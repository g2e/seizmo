function [data,stf]=reflect2seizmo(infile,outfile)
%REFLECT2SEIZMO    Constructs a SEIZMO dataset from Reflectivity I/O files
%
%    Usage:    data=reflect2seizmo(infile,outfile)
%
%    Description:
%     DATA=REFLECT2SEIZMO(INFILE,OUTFILE) converts reflectivity output to a
%     SEIZMO dataset using the input and output files from a reflect run.
%     This is poorly tested at the moment.
%
%    Notes:
%     - Why so slow?  Parsing ascii files in Matlab...
%     - Exponential factor:  data=data/e^(expfac*t)
%
%    Header changes: ALL
%
%    Examples:
%
%    See also: READ_REFLECT_INPUT, MAKE_REFLECT_INPUT, REFLECT2SEIZMO

%     Version History:
%        Aug. 10, 2010 - initial version
%        Feb.  2, 2011 - spread model name across KUSER0-2 to allow longer
%                        model names (one reason a fixed format fails)
%        Jan. 28, 2012 - minor improvement to strnlen usage
%        Jan. 26, 2014 - abs path exist fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2014 at 23:00 GMT

% todo:

% check nargin
error(nargchk(0,2,nargin));

% defaults
if(nargin<1); infile=[]; outfile=[]; end
if(nargin<2); outfile=[]; end

% read info
s=read_reflect_input(infile);

% attempt to read as binary if that fails then ascii
[x,p,raw]=read_reflect_binary_output(outfile,s.nsta,s.npts);

% this is for reading ascii output
%raw=read_reflect_output(outfile);
% parse the raw data
%raw=reshape(raw,[],s.nsta);
% first two points are distance & reduction slowness
%x=raw(1,:)'; p=raw(2,:)'; raw(1:2,:)=[];
% ZRTZRT...ZRT
%raw=reshape(raw,[],3*s.nsta);

% check distance
%[x s.dist x-s.dist abs(x-s.dist)>eps(single(s.dist))]
if(any(abs(x-s.dist)>100*eps(single(s.dist))))
    error('seizmo:reflect2seizmo:ioMismatch',...
        'Distance disagreement between input/output files!');
end

% find the CMT if available
% - this will provide the event info!
cmt=findcmt('name',s.event);
if(~isscalar(cmt.name))
    error('seizmo:reflect2seizmo:nonUniqueCMT',...
        'EVENT name does not resolve to a unique CMT!');
end

% get moment magnitude
mw=(2/3)*(log10(cmt.scalarmoment*10^cmt.exponent)-16.1);

% it is easy to deal with exponential factor right now
% - but I do not know if I'm doing this correctly (need an expert here)
expfac=exp((0:(s.npts-1))*s.delta*s.expfac)';
for i=1:3*s.nsta
    raw(:,i)=raw(:,i).*expfac;
end

% applying scale factors (not sure about this currently)
% this is the scale factor in reftosac (no explanation)
raw=-0.77e3*raw;
% this is more like what i expected
% cm ==> nm  = *1e4
% 1e25 ==> 1e27 = *1e2
raw=raw*1e4*10^(cmt.exponent-25);

% arrange into t,x pairings (blank t to conserve memory)
raw=[cell(1,3*s.nsta); mat2cell(raw,s.npts,ones(1,3*s.nsta))];

% now create seizmo records
% - we could do 3cmp records if I had a standard
data=bseizmo(raw{:});

% get some necessary values
idep='idisp'; if(s.units); idep='ivel'; end
time=fixtimes([cmt.year cmt.month cmt.day cmt.hour ...
    cmt.minute cmt.seconds+cmt.centroidtime],'utc');
kname=reshape(getwords(joinwords(s.staname,'.'),'.'),[],s.nsta)';
knetwk=kname(:,1); kstnm=kname(:,2);
khole=kname(:,3); kcmpnm=kname(:,4);
kcmpnm=strtrim(cellstr(strnlen(char(kcmpnm),2)));
[dist,az,baz]=sphericalinv(cmt.centroidlat,cmt.centroidlon,s.stla,s.stlo);
dist=dist*6371*pi/180;
b=s.vrstart+dist.*p;
e=b+(s.npts-1)*s.delta;

% check delaz
%[dist s.dist dist-s.dist abs(dist-s.dist)>eps(single(dist))]
if(any(abs(dist-s.dist)>100*eps(single(dist))))
    error('seizmo:reflect2seizmo:ioMismatch',...
        'Distance disagreement between input/cmt files!');
end
%[az s.az az-s.az abs(az-s.az)>eps(single(az))]
if(any(abs(az-s.az)>100*eps(single(az))))
    error('seizmo:reflect2seizmo:ioMismatch',...
        'Azimuth disagreement between input/cmt files!');
end

% split model string
s.model=strnlen(s.model,24);

% populate the headers
data(1:3:end)=changeheader(data(1:3:end),'b',b,'e',e,...
    'stla',s.stla,'stlo',s.stlo,'stel',0,'stdp',0,...
    'cmpaz',0,'cmpinc',0,...
    'knetwk',knetwk,'kstnm',kstnm,...
    'khole',khole,'kcmpnm',strcat(kcmpnm,'Z'));
data(2:3:end)=changeheader(data(2:3:end),'b',b,'e',e,...
    'stla',s.stla,'stlo',s.stlo,...
    'cmpaz',mod(baz-180,360),'cmpinc',90,...
    'knetwk',knetwk,'kstnm',kstnm,...
    'khole',khole,'kcmpnm',strcat(kcmpnm,'R'));
data(3:3:end)=changeheader(data(3:3:end),'b',b,'e',e,...
    'stla',s.stla,'stlo',s.stlo,...
    'cmpaz',mod(baz-90,360),'cmpinc',90,...
    'knetwk',knetwk,'kstnm',kstnm,...
    'khole',khole,'kcmpnm',strcat(kcmpnm,'T'));
data=changeheader(data,...
    'delta',s.delta,'isynth','ireflect',...
    'kevnm',s.event,'ievtyp','iquake','imagsrc','igcmt','ko','CMT',...
    'idep',idep,'iztype','io','o',0,'mag',mw,'imagtyp','imw','evel',0,...
    'evla',cmt.centroidlat,'evlo',cmt.centroidlon,...
    'evdp',cmt.centroiddep*1000,'stel',0,'stdp',0,'z6',{time},...
    'lcalda',true,'lpspol',true,'leven',true,'lovrok',true,...
    'kdatrd',s.run,'user0',s.nslow,'kuser0',s.model(1:8),...
    'user1',s.expfac,'kuser1',s.model(9:16),'kuser2',s.model(17:24),...
    'user3',s.freqlimits(1),'user4',s.freqlimits(2),...
    'user5',s.slowlimits(1),'user6',s.slowlimits(2),...
    'user7',s.taperlimits(1),'user8',s.taperlimits(2),...
    'resp0',s.filter,'resp1',s.field,'resp2',s.response,...
    'resp3',s.space,'resp4',s.function,'resp5',s.radiation,...
    'resp6',s.vr,'resp7',s.vrstart,...
    'resp8',s.dslimits(1),'resp9',s.dslimits(2));

% set name
data=genname(data);

% force check header
old=checkheader_state(true);
data=checkheader(data);
checkheader_state(old);

% source-time function
stf=s.srcfun;

end


function [x,p,raw]=read_reflect_binary_output(file,nsta,npts)
%READ_REFLECT_BINARY_OUTPUT  Reads in reflectivity binary output

% directory separator
fs=filesep;

% graphical selection
if(isempty(file))
    [file,path]=uigetfile(...
        {'*.sei;*.SEI' 'Reflect Output Files (*.sei,*.SEI)';
        '*.*' 'All Files (*.*)'},'Select File');
    if(isequal(0,file))
        error('seizmo:read_reflect_binary_output:noFileSelected',...
            'No input file selected!');
    end
    file=[path fs file];
else
    % check file
    if(~isstring(file))
        error('seizmo:read_reflect_binary_output:fileNotString',...
            'FILE must be a string!');
    end
    if(~isabspath(file)); file=[pwd fs file]; end
    if(~exist(file,'file'))
        error('seizmo:read_reflect_binary_output:fileDoesNotExist',...
            'File: %s\nDoes Not Exist!',file);
    elseif(exist(file,'dir'))
        error('seizmo:read_reflect_binary_output:dirConflict',...
            'File: %s\nIs A Directory!',file);
    end
end

% open file for reading
fid=fopen(file,'rb');

% check if file is openable
if(fid<0)
    error('seizmo:read_reflect_binary_output:cannotOpenFile',...
        'File: %s\nNot Openable!',file);
end

% loop over each station
raw=nan(3*npts+2,nsta);
for i=1:nsta
    fseek(fid,8,'cof'); % record start tag
    raw(:,i)=fread(fid,npts*3+2,'float32'); % x,p,data
    fseek(fid,8,'cof'); % record finish tag
end

% close file
fclose(fid);

% extract x,p
x=raw(1,:)';
p=raw(2,:)';
raw(1:2,:)=[];

% separate components: ZRTZRTZRT...
raw=reshape(raw,[],3*nsta);

end

