function [cmt]=findcmt(varargin)
%FINDCMT    Returns Global CMT that is closest to input time & position
%
%    Usage:    cmt=findcmt('option1',value1,...,'optionN',valueN)
%
%    Description: CMT=FINDCMT('OPTION1',VALUE1,...,'OPTIONN',VALUEN) finds
%     the Global CMT(s) (www.globalcmt.org) using the parameters provided.
%     The options are:
%      TIME      -- look for cmts near this time (see Notes for format)
%      LOCATION  -- look for cmts near this location ([deg_lat deg_lon])
%      DEPTH     -- look for cmts near this depth (in kilometers)
%      MAGNITUDE -- look for cmts near this moment magnitude
%      TIMESCALE -- scales time discrepancy to distance (default is 5 km/s)
%      MAGSCALE  -- scales magnitude discrepancy to distance (see Notes)
%      CATALOG   -- 'both' (default),'full', 'quick'
%      NUMCMT    -- integer>0 (default is 1)
%
%    Notes:
%     - TIME formats: [YEAR JULDAY] (start of day)
%                     [YEAR MONTH CALDAY] (start of day)
%                     [YEAR JULDAY HOUR MINUTE SECONDS]
%                     [YEAR MONTH CALDAY HOUR MINUTE SECONDS]
%     - MAGSCALE: 10^(MAGSCALE*abs(Mw-MAGNITUDE)) km (default is 4)
%     - CATALOG may be an NDK struct
%     - GlobalCMT catalogs are cached under global SEIZMO.GLOBALCMT
%
%    Examples:
%     Find the event closest to the new millinium:
%      findcmt('time',[2000 1])
%
%     Find deepest cmt:
%      findcmt('depth',1000)
%
%     Find cmt nearest the south pole:
%      findcmt('location',[-90 0])
%
%     Find largest cmt:
%      findcmt('magnitude',10)
%
%    See also: SETEVENT, READNDK, SSIDX

%     Version History:
%        Mar. 30, 2010 - initial version
%        July 30, 2010 - update for new ndk structs
%        Aug.  2, 2010 - catalog caching, both catalogs option
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  2, 2010 at 21:30 GMT

% todo:

% check optional inputs
opt=check_findcmt_options(varargin{:});

% get distances
dd=0;
if(~isempty(opt.DEPTH))
    dd=abs(opt.CATALOG.centroiddep-opt.DEPTH);
end
ld=0;
if(~isempty(opt.LOCATION))
    ld=vincentyinv(opt.LOCATION(1),opt.LOCATION(2),...
        opt.CATALOG.centroidlat,opt.CATALOG.centroidlon);
end
td=0;
if(~isempty(opt.TIME))
    td=opt.TIMESCALE*abs(timediff(opt.TIME,...
        [opt.CATALOG.year opt.CATALOG.month opt.CATALOG.day ...
        opt.CATALOG.hour opt.CATALOG.minute ...
        opt.CATALOG.seconds+opt.CATALOG.centroidtime]));
end
md=0;
if(~isempty(opt.MAGNITUDE))
    md=10.^(opt.MAGSCALE*...
        abs((2/3)*log10(opt.CATALOG.scalarmoment...
        .*10.^opt.CATALOG.exponent)-10.7-opt.MAGNITUDE));
end

% sort by cumulative distance and pull asked amount
d=sqrt(dd.^2+ld.^2+td.^2+md.^2);
[d,idx]=sort(d);
cmt=ssidx(opt.CATALOG,idx(1:opt.NUMCMT));

end


function [opt]=check_findcmt_options(varargin)

% global SEIZMO access for catalog caching
global SEIZMO

% require option/value pairings
if(mod(nargin,2))
    error('seizmo:findcmt:badInput',...
        'Unpaired OPTION/VALUE!');
end

% require options are strings
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:findcmt:badInput',...
        'OPTION must be a string!');
end

% defaults (not checked, so must be good)
opt.TIME=[]; % optional
opt.LOCATION=[]; % optional
opt.DEPTH=[]; % optional
opt.MAGNITUDE=[]; % optional
opt.TIMESCALE=5; % x km/s
opt.MAGSCALE=4; % 10^(y*abs(mw-mag))
opt.CATALOG='full'; % full/quick
opt.NUMCMT=1; % number of cmts returned

% valid strings
valid.CATALOG={'full' 'quick' 'both'};
valid.TIMES={'[YEAR JDAY]' '[YEAR MONTH CDAY]' ...
    '[YEAR JDAY HOUR MINUTE SECONDS]' ...
    '[YEAR MONTH CDAY HOUR MINUTE SECONDS]'};
valid.NDKFIELDS={'scalarmoment' 'exponent' 'year' 'month' 'day' 'hour' ...
    'minute' 'seconds' 'centroidtime' 'centroidlat' 'centroidlon' ...
    'centroiddep'};

% loop over options
for i=1:2:nargin
    switch lower(varargin{i})
        case {'time' 't'}
            % check
            sz=size(varargin{i+1});
            if(~isreal(varargin{i+1}) || numel(sz)>2 ...
                    || sz(1)~=1 || ~any(sz(2)==[2 3 5 6]))
                error('seizmo:findcmt:badInput',...
                    ['TIMES must have one of the ' ...
                    'following datetime formats:\n' ...
                    sprintf('%s\n',valid.TIMES)]);
            end
            opt.TIME=varargin{i+1};
        case {'location' 'position' 'loc' 'pos' 'l'}
            % check
            if(~isreal(varargin{i+1}) ...
                    || ~isequal(size(varargin{i+1}),[1 2]))
                error('seizmo:findcmt:badInput',...
                    'LOCATION must be in the format [LAT LON]!');
            end
            opt.LOCATION=varargin{i+1};
        case {'depth' 'dep' 'd' 'z'}
            % check
            if(~isreal(varargin{i+1}) || ~isscalar(varargin{i+1}))
                error('seizmo:findcmt:badInput',...
                    'DEPTH must be a real-valued scalar in kilometers!');
            end
            % auto-fix meter depths
            if(varargin{i+1}>1000)
                warning('seizmo:findcmt:badDepth',...
                    ['Assuming DEPTH>1000 is in meters.\n' ...
                    'Converting to kilometers!']);
                varargin{i+1}=varargin{i+1}/1000;
            end
            opt.DEPTH=varargin{i+1};
        case {'catalog' 'cat' 'c'}
            % check
            if(isstruct(varargin{i+1}) && all(ismember(valid.NDKFIELDS,...
                    fieldnames(varargin{i+1}))))
                opt.CATALOG=varargin{i+1};
            elseif(ischar(varargin{i+1}) ...
                    && ismember(lower(varargin{i+1}),valid.CATALOG))
                opt.CATALOG=varargin{i+1};
            else
                error('seizmo:findcmt:badInput',...
                    ['CATALOG must be a ndk-struct or:\n' ...
                    sprintf('''%s'' ',valid.CATALOG{:})]);
            end
        case {'timescale' 'ts' 'vel'}
            % check
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('seizmo:findcmt:badInput',...
                    'TIMESCALE must be a real-valued scalar in km/s!');
            end
            opt.TIMESCALE=varargin{i+1};
        case {'magscale' 'ms'}
            % check
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('seizmo:findcmt:badInput',...
                    'MAGSCALE must be a real-valued scalar!');
            end
            opt.MAGSCALE=varargin{i+1};
        case {'magnitude' 'mag' 'm'}
            % check
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('seizmo:findcmt:badInput',...
                    'MAGNITUDE must be a real-valued scalar!');
            end
            opt.MAGNITUDE=varargin{i+1};
        case {'numcmt' 'num' 'n'}
            % check
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}) ...
                    || fix(varargin{i+1})~=varargin{i+1} ...
                    || varargin{i+1}<1)
                error('seizmo:findcmt:badInput',...
                    'NUMCMT must be an integer >0!');
            end
            opt.NUMCMT=varargin{i+1};
        otherwise
            error('seizmo:findcmt:badInput',...
                'Unknown Option: %s !',varargin{i});
    end
end

% get catalog
if(ischar(opt.CATALOG))
    switch lower(opt.CATALOG)
        case 'both'
            % try loading cached version
            if(seizmoverbose)
                disp('Loading Catalog (May take a second...)');
            end
            try
                tmp1=SEIZMO.GLOBALCMT.FULL;
            catch
                % not there so load and cache
                try
                    tmp1=load('globalcmt_full');
                    SEIZMO.GLOBALCMT.FULL=tmp1;
                catch
                    error('seizmo:findcmt:badCatalog',...
                        'GlobalCMT catalogs are misconfigured!');
                end
            end
            try
                tmp2=SEIZMO.GLOBALCMT.QUICK;
            catch
                % not there so load and cache
                try
                    tmp2=load('globalcmt_quick');
                    SEIZMO.GLOBALCMT.QUICK=tmp2;
                catch
                    error('seizmo:findcmt:badCatalog',...
                        'GlobalCMT catalogs are misconfigured!');
                end
            end
            if(isstruct(tmp1) && ~all(ismember(valid.NDKFIELDS,...
                    fieldnames(tmp1))))
                error('seizmo:findcmt:badCatalog',...
                    'GlobalCMT catalog is corrupt!');
            end
            if(isstruct(tmp2) && ~all(ismember(valid.NDKFIELDS,...
                    fieldnames(tmp2))))
                error('seizmo:findcmt:badCatalog',...
                    'GlobalCMT catalog is corrupt!');
            end
            fields=fieldnames(tmp1);
            if(~isequal(fields,fieldnames(tmp2)))
                error('seizmo:findcmt:badCatalog',...
                    'GlobalCMT catalogs are misconfigured!');
            end
            opt.CATALOG=[];
            for i=1:numel(fields)
                opt.CATALOG.(fields{i})=...
                    [tmp1.(fields{i}); tmp2.(fields{i})];
            end
            if(seizmoverbose); disp('Loaded Catalog!'); end
        case {'full' 'quick'}
            % try loading cached version
            try
                opt.CATALOG=SEIZMO.GLOBALCMT.(upper(opt.CATALOG));
            catch
                % not there so load and cache
                try
                    if(seizmoverbose)
                        disp('Loading Catalog (May take a second...)');
                    end
                    cattype=upper(opt.CATALOG);
                    opt.CATALOG=load(['globalcmt_' opt.CATALOG]);
                    SEIZMO.GLOBALCMT.(cattype)=opt.CATALOG;
                    if(seizmoverbose); disp('Loaded Catalog!'); end
                catch
                    error('seizmo:findcmt:badCatalog',...
                        'GlobalCMT catalogs are misconfigured!');
                end
            end
            if(isstruct(opt.CATALOG) && ~all(ismember(valid.NDKFIELDS,...
                    fieldnames(opt.CATALOG))))
                error('seizmo:findcmt:badCatalog',...
                    'GlobalCMT catalog is corrupt!');
            end
    end
end

end

