function [varargout]=make_reflect_input(...
    run,cmt,model,kname,stla,stlo,varargin)
%MAKE_REFLECT_INPUT    Makes input file for a reflectivity run
%
%    Usage:    make_reflect_input(run,cmt,1dmodel,kname,stla,stlo)
%              make_reflect_input(...,'option1',val1,...'optionN',valN,...)
%              o=make_reflect_input(...)
%
%    Description:
%     MAKE_REFLECT_INPUT(RUN,CMT,1DMODEL,KNAME,STLA,STLO) will create
%     an input file for generating reflectivity synthetics with Brian
%     Kennett's Fortran-based reflect codes.  See the Notes section on
%     where to acquire the code that can run these input files (they are
%     non-standard -- but who's reflect codes are? :) ).  RUN is the a name
%     for this reflectivity run (like 'test001' etc).  CMT is a cmt struct
%     as output by FINDCMT and defines the event.  1DMODEL is the 1D
%     velocity model as output from for example PREM.  Note that AK135 &
%     IASP91 do not have a Q model included and so will need to be
%     augmented.  The KNAME, STLA, STLO arguments define the stations and
%     are all SEIZMO header fields.  KNAME is a Nx4 cell array of network,
%     station, location, and component codes.  A menu is presented to
%     indicate the output ascii file.
%
%     MAKE_REFLECT_INPUT(...,'OPTION1',VAL1,...'OPTIONN',VALN,...) sets
%     specific reflectivity options using the given values.  Available
%     options are as follows:
%      OPTION          VALID VALUES
%       file            valid file string or [] (gui selection)
%       field           'both' 'farfield' 'exactff'
%       response        'full' 'p' 's'
%       space           'freehalf' 'fohalf' 'freewhole' 'whole'
%       debug           true/false
%       function        'bessel' 'hankel'
%       radiation       'full' 'down'
%       srctype         'delta' 'gaussian' 'triangle' (sets custom to [])
%       customsrcfun    numeric vector
%       reflections     'none' 'primary' 'some' 'all'
%       delta           positive real
%       freqlimits      [lo hi] in +Hz
%       npts            positive integer corresponding to power of 2
%       expfac          real number
%       filter          'none' 'lp' 'bp' 'sinc'
%       units           'disp' 'velo'
%       nslow           positive integer
%       slowlimits      [lo hi] in s/km
%       denseslowlimits [lo hi] in s/km
%       taperlimits     [lo hi] in s/km
%       vr              real number in km/s
%       vrstart         real number in seconds
%
%      OPTION          DEFAULT VALUE
%       file            [] (gui select)
%       field           'exactff'
%       response        'full'
%       space           'freehalf'
%       debug           false
%       function        'bessel'
%       radiation       'full'
%       srctype         'delta'
%       customsrcfun    []
%       reflections     'all'
%       delta           1 second
%       freqlimits      [0 .01] Hz
%       npts            8192
%       expfac          0
%       filter          'none'
%       units           'disp'
%       nslow           4000
%       slowlimits      [0 .5] s/km
%       denseslowlimits []
%       taperlimits     []
%       vr              0 km/s
%       vrstart         0 seconds
%
%     O=MAKE_REFLECT_INPUT(...) returns an options struct useful for
%     debugging purposes.
%
%    Notes:
%     - Reflect code can be found in the same directory as this mfile in
%       the gzipped tarball reflect.tar.gz.  Please don't expect me to get
%       this running for you.
%
%    Examples:
%     % generate a reflectivity input file for a dataset
%     cmt=findcmt;
%     model=prem('r',[0 5000]); % reflectivity code doesn't handle icore
%     [kname,stla,stlo]=getheader(data,'kname','stla','stlo');
%     make_reflect_input('test',cmt,model,kname,stla,stlo);
%
%    See also: READ_REFLECT_INPUT, READ_REFLECT_OUTPUT, REFLECT2SEIZMO

%     Version History:
%        Aug. 10, 2010 - initial version
%        Feb. 28, 2011 - fix overwrite crash
%        Jan. 26, 2014 - abs path exist fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2014 at 23:00 GMT

% todo:

% check nargin
error(nargchk(6,inf,nargin));

% simple constants
Re=6371;
D2R=pi/180;

% check runname
if(~ischar(run) || isempty(run) || ndims(run)~=2 ...
        || size(run,1)~=1 || numel(getwords(run))~=1)
    error('seizmo:make_reflect_input:badInput',...
        'RUNNAME must be a string (and one word)!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% begin ndk section!!! %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% require certain fields exist
reqfields={'centroidlat' 'centroidlon' 'centroiddep' 'srcfuncdur' ...
    'mrr' 'mtt' 'mpp' 'mrt' 'mrp' 'mtp'};
if(any(~ismember(reqfields,fieldnames(cmt))))
    error('seizmo:make_reflect_input:badInput',...
        'CMT struct is insufficient!');
end

% require fields are scalar
for i=1:numel(reqfields)
    if(~isscalar(cmt.(reqfields{i})))
        error('seizmo:make_reflect_input:badInput',...
            'CMT struct field %s must be scalar!',reqfields{i});
    end
end

% get flattened source depth
fcmtdepth=flattendepth(cmt.centroiddep);

% get moment tensor
mt=[cmt.mrr cmt.mtt cmt.mpp cmt.mrt cmt.mrp cmt.mtp];
mt(5:6)=-mt(5:6);
mt=[mt([2 6 4]); mt([6 3 5]); mt([4 5 1])];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% end ndk section!!! %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% begin model section!!! %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check model
error(chk1dmodel(model));

% require at least 2 depths
if(numel(model.depth)<2)
    error('seizmo:make_reflect_input:badInput',...
        '1DMODEL must have 2 or more depth points!');
end

% must be isotropic
if(~model.isotropic)
    error('seizmo:make_reflect_input:badInput',...
        '1DMODEL must be isotropic!');
end

% require vp, vs, rho
fields=fieldnames(model);
ok=all(ismember({'vp' 'vs' 'rho'},fields));
if(~ok)
    error('seizmo:make_reflect_input:badInput',...
        '1DMODEL must have vp, vs & rho properties!');
end

% check for solid below liquid
if(any(diff(model.vs==0)==-1))
    warning('seizmo:make_reflect_input:badInput',...
        '1DMODEL has solid layers below liquid layers!');
end

% require model to have Q fields
qku=all(ismember({'qk' 'qu'},fields));
qps=all(ismember({'qp' 'qs'},fields));
if(~qku && ~qps)
    error('seizmo:make_reflect_input:badInput',...
        '1DMODEL must have Q properties!');
elseif(qku && ~qps)
    % convert to qp & qs
    model.qp=qkqu2qp(model.qk,model.qu,model.vp,model.vs);
    model.qs=model.qu;
end

% flatten if necessary
if(~model.flattened)
    model=flatten_1dmodel(model);
end

% tabularize model (aka make layers of uniform properties)
model=tabularize_1dmodel(model);

% convert qp & qs to attenuation
model.qp=1./model.qp;
model.qs=1./model.qs;

% number of model layers
nlayers=numel(model.thick);

% what layer is the cmt in?
cmtlayer=find(fcmtdepth>model.depth,1,'last');
if(isempty(cmtlayer) || cmtlayer>nlayers)
    error('seizmo:make_reflect_input:badInput',...
        'CMT depth is outside model tabulation!');
end

% make model matrix
mm=[model.vp model.vs model.rho model.thick model.qp model.qs];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% end model section!!! %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% begin station section!!! %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check kname
if(~iscellstr(kname) || ndims(kname)~=2 || size(kname,2)~=4)
    error('seizmo:make_reflect_input:badInput',...
        'KNAME must be an Nx4 cellstr array!');
end

% combine kname
kname=strcat(kname(:,1),'.',kname(:,2),'.',kname(:,3),'.',kname(:,4));

% check stla/stlo
if(~isreal(stla) || ~isreal(stlo))
    error('seizmo:make_reflect_input:badInput',...
        'STLA & STLO must be real valued arrays!');
elseif(~isequalsizeorscalar(kname,stla(:),stlo(:)))
    error('seizmo:make_reflect_input:badInput',...
        'KNAME, STLA, STLO must be equal length or scalar!');
end
[kname,stla,stlo]=expandscalars(kname,stla,stlo);
nsta=numel(kname);

% fix stla/stlo
[stla,stlo]=fixlatlon(stla,stlo);

% get dist/az from stla/stlo & cmt location
[dist,az]=sphericalinv(cmt.centroidlat,cmt.centroidlon,stla,stlo);
dist=dist*Re*D2R;

% set variable azimuth flag to true
vazflag=true;

% create stainfo matrix
% name dist azi lat lon
stainfo=[kname num2cell([dist az stla stlo])]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% end station section!!! %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% begin optional section!!! %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse optional inputs
o=parse_reflect_options(varargin{:});

% prefix reflection level to model matrix
mm=[o.reflections*ones(nlayers,1) mm];

% handle source time function
switch lower(o.srctype)
    case 'delta'
        srcfun=1;
    case 'custom'
        srcfun=o.custom;
    case {'gaussian' 'triangle'}
        srcfun=make_source_timefunction(o.delta,cmt.srcfuncdur,o.srctype);
        srcfun=srcfun{1};
    otherwise
end
nsrcfun=numel(srcfun);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% end optional section!!! %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% begin write section!!! %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% directory separator
fs=filesep;

% get output filename if one not given
if(isempty(o.file))
    [o.file,path]=uiputfile(...
        {'*.DAT;*.dat' 'Reflect Input Files (*.DAT,*.dat)';
        '*.*' 'All Files (*.*)'},...
        'Save reflect input as',[run '.dat']);
    if(isequal(0,o.file))
        error('seizmo:make_reflect_input:noFileSelected',...
            'No output file selected!');
    end
    o.file=[path fs o.file];
else
    % check file
    if(~isstring(o.file))
        error('seizmo:make_reflect_input:fileNotString',...
            'FILE must be a string!');
    end
    if(~isabspath(o.file)); o.file=[pwd fs o.file]; end
    if(exist(o.file,'file'))
        if(exist(o.file,'dir'))
            error('seizmo:make_reflect_input:dirConflict',...
                'Reflect Input File: %s\nIs A Directory!',o.file);
        end
        if(false)
            fprintf('Reflect Input File: %s\nFile Exists!\n',o.file);
            reply=input('Overwrite? Y/N [N]: ','s');
            if(isempty(reply) || ~strncmpi(reply,'y',1))
                disp('Not overwriting!');
                return;
            end
            disp('Overwriting!');
        end
    end
end

% open file for writing as ascii
fid=fopen(o.file,'wt');

% check if file is openable
if(fid<0)
    error('seizmo:make_reflect_input:cannotOpenFile',...
        'Reflect Input File: %s\nNot Openable!',file);
end

% print to file
fprintf(fid,'%-80s\n',run);
fprintf(fid,'%-80s\n',cmt.name{1});
fprintf(fid,'%-80s\n',model.name);
fprintf(fid,'%4i',nlayers,nsrcfun,nsta,o.field,o.response,o.space,...
    o.debug,o.function,o.radiation); fprintf(fid,'\n');
fprintf(fid,'%10i%10i%10.5f%10.5f\n',o.npts,log2(o.npts),o.delta,o.expfac);
fprintf(fid,'%10i%10.5f\n',cmtlayer,fcmtdepth);
fprintf(fid,'%10.5f%10.5f%10.5f\n',mt');
fprintf(fid,'%10.5f%10.5f\n',o.freqlimits);
fprintf(fid,'%10i%10.5f%10.5f\n',o.nslow,o.slowlimits);
fprintf(fid,'%10i%10.5f%10.5f\n',o.dsflag,o.dslimits);
fprintf(fid,'%10i%10.5f%10.5f\n',o.taperflag,o.taperlimits);
fprintf(fid,'%i%9.2f%10.5f%10.5f%10.5f%10.5f%10.5f\n',mm');
fprintf(fid,'%4i%4i%4i\n',vazflag,o.filter,o.units);
fprintf(fid,'%-20s%10.3f%10.5f%10.5f%10.5f\n',stainfo{:});
fprintf(fid,'%10.5f',[real(srcfun); imag(srcfun)]); fprintf(fid,'\n');
fprintf(fid,'%10.5f%10.5f\n',o.vr,o.vrstart);

% close file
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% end write section!!! %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargout)
    varargout{1}=o;
end

end


function [model]=tabularize_1dmodel(model)
% is bottom layer a discontinuity? if so then repeat bottom point
dcbot=model.depth(end)-model.depth(end-1)==0;
if(dcbot)
    model.depth=model.depth([1:end end]);
    model.vp=model.vp([1:end end]);
    model.vs=model.vs([1:end end]);
    model.rho=model.rho([1:end end]);
    model.qp=model.qp([1:end end]);
    model.qs=model.qs([1:end end]);
end
% thicknesses & averages
model.thick=model.depth(2:end)-model.depth(1:end-1);
model.vp=(model.vp(2:end)+model.vp(1:end-1))./2;
model.vs=(model.vs(2:end)+model.vs(1:end-1))./2;
model.rho=(model.rho(2:end)+model.rho(1:end-1))./2;
model.qp=(model.qp(2:end)+model.qp(1:end-1))./2;
model.qs=(model.qs(2:end)+model.qs(1:end-1))./2;
% toss discontinuity layers (except last)
dcidx=find(model.thick(1:end-1)==0);
model.thick(dcidx)=[];
model.vp(dcidx)=[];
model.vs(dcidx)=[];
model.rho(dcidx)=[];
model.qp(dcidx)=[];
model.qs(dcidx)=[];
% repeat top layer (but with 0 thickness)
model.thick=[0; model.thick];
model.vp=model.vp([1 1:end]);
model.vs=model.vs([1 1:end]);
model.rho=model.rho([1 1:end]);
model.qp=model.qp([1 1:end]);
model.qs=model.qs([1 1:end]);
end


function [depth]=flattendepth(depth)
Re=6371;
depth=-Re*log((Re-depth)/Re);
end


%function [depth]=unflattendepth(depth)
%Re=6371;
%depth=Re*(1-exp(-depth/Re));
%end


function [o]=parse_reflect_options(varargin)

% set defaults
varargin=[{'field' 'exactff' 'resp' 'full' 'space' 'halffree' 'dt' 1 ...
    'bug' false 'func' 'bessel' 'rad' 'full' 'src' 'delta' 'refl' 'all' ...
    'npts' 8192 'freq' [0 0.01] 'ef' 0 'filt' 'none' 'u' 'disp' ...
    'ns' 4000 'sl' [0 0.5] 'dsl' [] 'tl' [] 'vr' 0 'vrs' 0 'file' []} ...
    varargin];

% require option/value pairs
if(mod(numel(varargin),2))
    error('seizmo:make_reflect_input:badInput',...
        'Unpaired Option/Value!');
end

% require options are strings
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:make_reflect_input:badInput',...
        'Options must be specified by a string!');
end

% loop over every pair
for i=1:2:numel(varargin)
    % parse optional inputs
    v=varargin{i+1};
    switch lower(varargin{i})
        case {'file'}
            % o.file
            if(isempty(v))
                o.file=[];
            elseif(~ischar(v) || ndims(v)~=2 || size(v,1)~=1)
                error('seizmo:make_reflect_input:badInput',...
                    'FILE option must be a string!');
            end
            o.file=v;
        case {'field' 'nff'}
            % o.field ('both' 'ff' 'exactff')
            if(isempty(v)); continue; end
            if(~ischar(v))
                error('seizmo:make_reflect_input:badInput',...
                    'FIELD option must be a string!');
            end
            switch lower(v)
                case {'b' 'both'}
                    o.field=0;
                case {'ff' 'farfield'}
                    o.field=1;
                case {'eff' 'exactff'}
                    o.field=2;
                otherwise
                    error('seizmo:make_reflect_input:badInput',...
                        'Unknown FIELD type: %s',v);
            end
        case {'response' 'resp'}
            % o.response ('full' 'p' 's')
            if(isempty(v)); continue; end
            if(~ischar(v))
                error('seizmo:make_reflect_input:badInput',...
                    'RESPONSE option must be a string!');
            end
            switch lower(v)
                case {'f' 'full'}
                    o.response=0;
                case {'p'}
                    o.response=1;
                case {'s'}
                    o.response=2;
                otherwise
                    error('seizmo:make_reflect_input:badInput',...
                        'Unknown RESPONSE type: %s',v);
            end
        case {'space' 'sp'}
            % o.space ('halffree' 'halffo' 'wholefree' 'whole')
            if(isempty(v)); continue; end
            if(~ischar(v))
                error('seizmo:make_reflect_input:badInput',...
                    'SPACE option must be a string!');
            end
            switch lower(v)
                case {'hf' 'fh' 'halffree' 'freehalf' 'freehalfspace'}
                    o.space=0;
                case {'hfo' 'foh' 'fohalf' 'halffo' 'firstorderhalfspace'}
                    o.space=1;
                case {'wf' 'fw' 'freewhole' 'wholefree' 'freewholespace'}
                    o.space=2;
                case {'whole' 'w'}
                    o.space=3;
                otherwise
                    error('seizmo:make_reflect_input:badInput',...
                        'Unknown SPACE type: %s',v);
            end
        case {'debug' 'bug' 'de'}
            % o.debug (true/false)
            if(isempty(v)); continue; end
            if(~islogical(v) || ~isreal(v) || ~isscalar(v))
                error('seizmo:make_reflect_input:badInput',...
                    'DEBUG option must be TRUE/FALSE!');
            end
            o.debug=v;
        case {'function' 'func'}
            % o.function ('bessel' 'hankel')
            if(isempty(v)); continue; end
            if(~ischar(v))
                error('seizmo:make_reflect_input:badInput',...
                    'FUNCTION option must be a string!');
            end
            switch lower(v)
                case {'b' 'bessel'}
                    o.function=0;
                case {'h' 'hankel'}
                    o.function=1;
                otherwise
                    error('seizmo:make_reflect_input:badInput',...
                        'Unknown FUNCTION type: %s',v);
            end
        case {'radiation' 'radpat' 'rad'}
            % o.radiation ('full' 'down')
            if(isempty(v)); continue; end
            if(~ischar(v))
                error('seizmo:make_reflect_input:badInput',...
                    'RADIATION option must be a string!');
            end
            switch lower(v)
                case {'f' 'full'}
                    o.radiation=0;
                case {'d' 'down'}
                    o.radiation=1;
                otherwise
                    error('seizmo:make_reflect_input:badInput',...
                        'Unknown RADIATION type: %s',v);
            end
        case {'custom' 'customsrcfun' 'cu' 'srcfun'}
            % o.custom (default empty ==> sets srctype to 'custom')
            if(isempty(v) || ~isnumeric(v) || ~isvector(v))
                error('seizmo:make_reflect_input:badInput',...
                    'CUSTOMSRCFUN must be a numeric vector!');
            else
                o.srctype='custom';
                o.custom=v(:)';
            end
        case {'sourcetype' 'srctype' 'src'}
            % o.srctype ('delta' 'gaussian' 'triangle') (sets custom to [])
            if(isempty(v)); continue; end
            if(~ischar(v))
                error('seizmo:make_reflect_input:badInput',...
                    'SOURCETYPE option must be a string!');
            end
            switch lower(v)
                case {'d' 'delta'}
                    o.srctype='delta';
                    o.custom=[];
                case {'g' 'gauss' 'gaussian'}
                    o.srctype='gaussian';
                    o.custom=[];
                case {'t' 'tri' 'triangle'}
                    o.srctype='triangle';
                    o.custom=[];
                otherwise
                    error('seizmo:make_reflect_input:badInput',...
                        'Unknown SOURCETYPE type: %s',v);
            end
        case {'reflections' 'refl'}
            % o.reflections ('none' 'primary' 'some' 'all')
            if(isempty(v)); continue; end
            if(~ischar(v))
                error('seizmo:make_reflect_input:badInput',...
                    'REFRECTIONS option must be a string!');
            end
            switch lower(v)
                case {'n' 'none'}
                    o.reflections=0;
                case {'p' 'primary'}
                    o.reflections=1;
                case {'s' 'some'}
                    o.reflections=2;
                case {'a' 'all'}
                    o.reflections=3;
                otherwise
                    error('seizmo:make_reflect_input:badInput',...
                        'Unknown REFLECTIONS type: %s',v);
            end
        case {'delta' 'dt'}
            % o.delta (sample spacing in seconds)
            if(isempty(v)); continue; end
            if(~isreal(v) || ~isscalar(v) || v<=0)
                error('seizmo:make_reflect_input:badInput',...
                    'DELTA option must be a positive scalar in seconds!');
            end
            o.delta=v;
        case {'npts' 'np'}
            % o.npts (mod(log2(npts),1) must be 0)
            if(isempty(v)); continue; end
            if(~isreal(v) || ~isscalar(v) || v~=fix(v) || v<1 ...
                    || mod(log2(v),1))
                error('seizmo:make_reflect_input:badInput',...
                    'NPTS must be a positive power of two (8192, etc)!');
            end
            o.npts=v;
        case {'freqlimits' 'freq' 'fl'}
            % o.freqlimits ([lo hi] in hz)
            if(isempty(v)); continue; end
            if(~isreal(v) || numel(v)~=2 || v(2)<=v(1) || any(v<0))
                error('seizmo:make_reflect_input:badInput',...
                    'FREQLIMITS must be [Lo Hi] in Hz!');
            end
            o.freqlimits=[v(1) v(2)];
        case {'expfac' 'exp' 'ef'}
            % o.expfac (must be zero if far field)
            if(isempty(v)); continue; end
            if(~isreal(v) || ~isscalar(v))
                error('seizmo:make_reflect_input:badInput',...
                    'EXPFAC option must be a real-valued scalar!');
            end
            o.expfac=v;
        case {'filter' 'filt'}
            % o.filter ('none' 'lp' 'bp' 'sinc')
            if(isempty(v)); continue; end
            if(~ischar(v))
                error('seizmo:make_reflect_input:badInput',...
                    'FILTER option must be a string!');
            end
            switch lower(v)
                case {'n' 'none'}
                    o.filter=0;
                case {'lp'}
                    o.filter=1;
                case {'bp'}
                    o.filter=2;
                case {'s' 'sinc'}
                    o.filter=4;
                otherwise
                    error('seizmo:make_reflect_input:badInput',...
                        'Unknown FILTER type: %s',v);
            end
        case {'units' 'u'}
            % o.units ('disp' 'velo')
            if(isempty(v)); continue; end
            if(~ischar(v))
                error('seizmo:make_reflect_input:badInput',...
                    'UNITS option must be a string!');
            end
            switch lower(v)
                case {'d' 'dis' 'disp'}
                    o.units=0;
                case {'v' 'vel' 'velo'}
                    o.units=1;
                otherwise
                    error('seizmo:make_reflect_input:badInput',...
                        'Unknown UNITS type: %s',v);
            end
        case {'nslow' 'ns'}
            % o.nslow (positive integer)
            if(isempty(v)); continue; end
            if(~isreal(v) || ~isscalar(v) || v~=fix(v) || v<1)
                error('seizmo:make_reflect_input:badInput',...
                    'NSLOW must be a positive integer!');
            end
            o.nslow=v;
        case {'slowlimits' 'sl'}
            % o.slowlimits ([lo hi] in s/km)
            if(isempty(v)); continue; end
            if(~isreal(v) || numel(v)~=2 || v(2)<=v(1) || any(v<0))
                error('seizmo:make_reflect_input:badInput',...
                    'SLOWLIMITS must be [Lo Hi] in s/km!');
            end
            o.slowlimits=[v(1) v(2)];
        case {'denseslowlimits' 'dslimits' 'dsl'}
            % o.denseslowlimits ([lo hi] in s/km) (sets dsflag)
            if(isempty(v))
                o.dsflag=false;
                o.dslimits=[0 0];
            elseif(~isreal(v) || numel(v)~=2 || v(2)<=v(1) || any(v<0))
                error('seizmo:make_reflect_input:badInput',...
                    'DENSESLOWLIMITS must be [Lo Hi] in s/km!');
            else
                o.dsflag=true;
                o.dslimits=[v(1) v(2)];
            end
        case {'slowtaper' 'taperlimits' 'taper' 'tl'}
            % o.taperlimits ([lo hi] in s/km) (sets slowtaper)
            if(isempty(v))
                o.taperflag=false;
                o.taperlimits=[0 0];
            elseif(~isreal(v) || numel(v)~=2 || v(2)<=v(1) || any(v<0))
                error('seizmo:make_reflect_input:badInput',...
                    'TAPERLIMITS must be [Lo Hi] in s/km!');
            else
                o.taperflag=true;
                o.taperlimits=[v(1) v(2)];
            end
        case {'vr' 'vred' 'redvelo'}
            % o.vr (km/s)
            if(isempty(v)); continue; end
            if(~isreal(v) || ~isscalar(v))
                error('seizmo:make_reflect_input:badInput',...
                    'VR option must be a scalar in km/s!');
            end
            o.vr=v;
        case {'vrstart' 'vrs' 'stmin'}
            % o.vrstart (seconds)
            if(isempty(v)); continue; end
            if(~isreal(v) || ~isscalar(v))
                error('seizmo:make_reflect_input:badInput',...
                    'VRSTART must be a real-valued scalar in seconds!');
            end
            o.vrstart=v;
        otherwise
            error('seizmo:make_reflect_input:badInput',...
                'Unknown Reflect Option: %s',varargin{i});
    end
end

% check that expfac is 0 if ff
if(o.expfac && o.field)
%    error('seizmo:make_reflect_input:badInput',...
%        'EXPFAC must be 0 if FIELD is of a far-field type!');
end

end
