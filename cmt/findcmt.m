function [cmt,d]=findcmt(varargin)
%FINDCMT    Returns Global CMT that is closest to input time & position
%
%    Usage:    cmt=findcmt('option1',value1,...,'optionN',valueN)
%              [cmt,dist]=findcmt(...)
%
%    Description:
%     CMT=FINDCMT('OPTION1',VALUE1,...,'OPTIONN',VALUEN) finds the Global
%     CMT(s) (www.globalcmt.org) using the parameters provided.  The
%     options are:
%      TYPE      -- location type: 'centroid' or 'hypocenter' (default)
%      TIME      -- look for cmts near this time (see Notes for format)
%      LOCATION  -- look for cmts near this location ([deg_lat deg_lon])
%      DEPTH     -- look for cmts near this depth (in kilometers)
%      MAGNITUDE -- look for cmts near this magnitude
%      MAGTYPE   -- magnitude type: 'mb' 'ms' or 'mw' (default)
%      TIMESCALE -- scales time discrepancy to distance (default is 5 km/s)
%      MAGSCALE  -- scales magnitude discrepancy to distance (see Notes)
%      CATALOG   -- CMT catalog: 'Full', 'Quick', 'Both' (Default)
%      NUMCMT    -- returns this many cmts, integer>0 (default is 1)
%      NAME      -- limits search to cmt(s) who's name satisfies the given
%                   regular expression (see REGEXP for details)
%
%     [CMT,DIST]=FINDCMT(...) returns the total distance factor used in
%     determining the best CMT for the inputs given.
%
%    Notes:
%     - Using the NAME option ignores all options but CATALOG
%     - TIME formats: [YEAR JULDAY] (start of day)
%                     [YEAR MONTH CALDAY] (start of day)
%                     [YEAR JULDAY HOUR MINUTE SECONDS]
%                     [YEAR MONTH CALDAY HOUR MINUTE SECONDS]
%     - MAGSCALE: MAGSCALE*abs(Mw-MAGNITUDE) km (default is 100)
%     - CATALOG may be an NDK struct
%     - GlobalCMT catalogs are cached under global SEIZMO.GLOBALCMT
%
%    Examples:
%     % Find the event closest to the new millinium:
%     findcmt('time',[2000 1])
%
%     % Find deepest cmt:
%     findcmt('depth',1000)
%
%     % Find cmt nearest the south pole:
%     findcmt('location',[-90 0])
%
%     % Find largest cmt:
%     findcmt('magnitude',10)
%
%     % Find largest cmt from Jan 2008:
%     findcmt('name','^\w200801','magnitude',10)
%
%    See also: FINDCMTS, SETEVENT, READNDK, SSIDX

%     Version History:
%        Mar. 30, 2010 - initial version
%        July 30, 2010 - update for new ndk structs
%        Aug.  2, 2010 - catalog caching, both catalogs option
%        Aug. 10, 2010 - add name option, handle large numcmt gracefully
%        Aug. 11, 2010 - altered/fixed magnitude scaling, added magtype opt
%        Sep. 28, 2010 - also output total distance
%        Apr. 28, 2011 - return n when n asked for and nothing else
%        Feb. 29, 2012 - doc update, type option & switch to hypo, use both
%                        catalogs by default
%        Mar.  1, 2012 - auto-creation of catalogs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2012 at 21:30 GMT

% todo:

% check optional inputs
opt=check_findcmt_options(varargin{:});

% name shortcut
if(~isempty(opt.NAME))
    idx=~cellfun('isempty',regexp(opt.CATALOG.name,opt.NAME));
    opt.CATALOG=ssidx(opt.CATALOG,idx);
    if(~sum(idx)); cmt=opt.CATALOG; d=[]; return; end
end

% type-dependant distances
switch lower(opt.TYPE)
    case {'centroid' 'cmt' 'c'}
        % get distances
        dd=0*opt.CATALOG.centroiddep;
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
    case {'hypocenter' 'hypo' 'h'}
        % get distances
        dd=0*opt.CATALOG.depth;
        if(~isempty(opt.DEPTH))
            dd=abs(opt.CATALOG.depth-opt.DEPTH);
        end
        ld=0;
        if(~isempty(opt.LOCATION))
            ld=vincentyinv(opt.LOCATION(1),opt.LOCATION(2),...
                opt.CATALOG.latitude,opt.CATALOG.longitude);
        end
        td=0;
        if(~isempty(opt.TIME))
            td=opt.TIMESCALE*abs(timediff(opt.TIME,...
                [opt.CATALOG.year opt.CATALOG.month opt.CATALOG.day ...
                opt.CATALOG.hour opt.CATALOG.minute opt.CATALOG.seconds]));
        end
end

% magnitude distance does not depend on type
md=0;
if(~isempty(opt.MAGNITUDE))
    switch lower(opt.MAGTYPE)
        case 'mw'
            m=(2/3).*(log10(opt.CATALOG.scalarmoment...
                .*10.^opt.CATALOG.exponent)-16.1);
        case 'ms'
            m=opt.CATALOG.ms;
        case 'mb'
            m=opt.CATALOG.ms;
    end
    md=opt.MAGSCALE.*abs(m-opt.MAGNITUDE);
end

% sort by cumulative distance and pull asked amount
d=sqrt(dd.^2+ld.^2+td.^2+md.^2);
[d,idx]=sort(d);
cmt=ssidx(opt.CATALOG,idx(1:min(numel(idx),opt.NUMCMT)));
d=d(1:min(numel(idx),opt.NUMCMT));

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
opt.NAME=[]; % optional
opt.TYPE='hypo'; % hypo, centroid
opt.TIME=[]; % optional
opt.LOCATION=[]; % optional
opt.DEPTH=[]; % optional
opt.MAGNITUDE=[]; % optional
opt.MAGTYPE='mw'; % ms, mb, mw
opt.TIMESCALE=5; % x km/s
opt.MAGSCALE=100; % km/mag
opt.CATALOG='both'; % full/quick/both
opt.NUMCMT=1; % number of cmts returned

% valid strings
valid.CATALOG={'full' 'quick' 'both'};
valid.TYPE={'hypocenter' 'hypo' 'h' 'centroid' 'cmt' 'c'};
valid.MAGTYPE={'mb' 'ms' 'mw'};
valid.TIMES={'[YEAR JDAY]' '[YEAR MONTH CDAY]' ...
    '[YEAR JDAY HOUR MINUTE SECONDS]' ...
    '[YEAR MONTH CDAY HOUR MINUTE SECONDS]'};
valid.NDKFIELDS={'scalarmoment' 'exponent' 'year' 'month' 'day' 'hour' ...
    'minute' 'seconds' 'centroidtime' 'centroidlat' 'centroidlon' ...
    'centroiddep' 'latitude' 'longitude' 'depth'};

% loop over options
for i=1:2:nargin
    switch lower(varargin{i})
        case {'name' 'na'}
            if(~ischar(varargin{i+1}) || ndims(varargin{i+1})~=2 ...
                    || size(varargin{i+1},1)~=1)
                error('seizmo:findcmt:badInput',...
                    'NAME must be a string!');
            end
            opt.NAME=varargin{i+1};
        case {'type' 'ty'}
            if(ischar(varargin{i+1}) ...
                    && ismember(lower(varargin{i+1}),valid.TYPE))
                opt.TYPE=varargin{i+1};
            else
                error('seizmo:findcmt:badInput',...
                    ['TYPE must be one of the following:\n' ...
                    sprintf('''%s'' ',valid.TYPE{:})]);
            end
        case {'time' 'ti' 't'}
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
        case {'magtype' 'mt'}
            if(ischar(varargin{i+1}) ...
                    && ismember(lower(varargin{i+1}),valid.MAGTYPE))
                opt.MAGTYPE=varargin{i+1};
            else
                error('seizmo:findcmt:badInput',...
                    ['MAGTYPE must be one of the following:\n' ...
                    sprintf('''%s'' ',valid.MAGTYPE{:})]);
            end
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
    opt.CATALOG=lower(opt.CATALOG);
    switch opt.CATALOG
        case 'both'
            % try loading cached version
            try
                tmp1=SEIZMO.GLOBALCMT.FULL;
            catch
                % no cache there...so load and cache
                try
                    if(seizmoverbose)
                        disp('Loading Catalog (May take a second...)');
                    end
                    tmp1=load('globalcmt_full');
                    SEIZMO.GLOBALCMT.FULL=tmp1;
                    if(seizmoverbose); disp('Loaded Catalog!'); end
                catch
                    % what happened?
                    if(exist('globalcmt_full.mat','file'))
                        % catalog is broken
                        disp('Local GlobalCMT Full catalog is corrupted!');
                        reply=input('Attempt to re-create? Y/N [N]: ','s');
                        if(~isempty(reply) && strncmpi(reply,'y',1))
                            globalcmt_create;
                            globalcmt_update;
                            tmp1=load('globalcmt_full');
                            SEIZMO.GLOBALCMT.FULL=tmp1;
                            if(seizmoverbose); disp('Loaded Catalog!'); end
                        else
                            error('seizmo:findcmt:badCatalog',...
                                'No Valid Local GlobalCMT Full catalog!');
                        end
                    else
                        % no catalog...ask to make one
                        disp('No Local GlobalCMT Full Catalog Found!');
                        reply=input('Attempt to create? Y/N [N]: ','s');
                        if(~isempty(reply) && strncmpi(reply,'y',1))
                            globalcmt_create;
                            globalcmt_update;
                            tmp1=load('globalcmt_full');
                            SEIZMO.GLOBALCMT.FULL=tmp1;
                            if(seizmoverbose); disp('Loaded Catalog!'); end
                        else
                            error('seizmo:findcmt:badCatalog',...
                                'No Local GlobalCMT Full catalog!');
                        end
                    end
                end
            end
            try
                tmp2=SEIZMO.GLOBALCMT.QUICK;
            catch
                % no cache there...so load and cache
                try
                    tmp2=load('globalcmt_quick');
                    SEIZMO.GLOBALCMT.QUICK=tmp2;
                catch
                    % what happened?
                    if(exist('globalcmt_quick.mat','file'))
                        % catalog is broken
                        disp('Local GlobalCMT Quick catalog corrupted!');
                        reply=input('Attempt to re-create? Y/N [N]: ','s');
                        if(~isempty(reply) && strncmpi(reply,'y',1))
                            globalcmt_update;
                            tmp2=load('globalcmt_quick');
                            SEIZMO.GLOBALCMT.QUICK=tmp2;
                            if(seizmoverbose); disp('Loaded Catalog!'); end
                        else
                            error('seizmo:findcmt:badCatalog',...
                                'No Valid Local GlobalCMT Quick catalog!');
                        end
                    else
                        % no catalog...ask to make one
                        disp('No Local GlobalCMT Quick Catalog Found!');
                        reply=input('Attempt to create? Y/N [N]: ','s');
                        if(~isempty(reply) && strncmpi(reply,'y',1))
                            globalcmt_update;
                            tmp2=load('globalcmt_quick');
                            SEIZMO.GLOBALCMT.QUICK=tmp2;
                            if(seizmoverbose); disp('Loaded Catalog!'); end
                        else
                            error('seizmo:findcmt:badCatalog',...
                                'No Local GlobalCMT Quick catalog!');
                        end
                    end
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
        case {'full' 'quick'}
            % try loading cached version
            try
                opt.CATALOG=SEIZMO.GLOBALCMT.(upper(opt.CATALOG));
            catch
                % no cache there...so load and cache
                try
                    if(seizmoverbose)
                        disp('Loading Catalog (May take a second...)');
                    end
                    cattype=upper(opt.CATALOG);
                    opt.CATALOG=load(['globalcmt_' opt.CATALOG]);
                    SEIZMO.GLOBALCMT.(cattype)=opt.CATALOG;
                    if(seizmoverbose); disp('Loaded Catalog!'); end
                catch
                    % what happened?
                    if(exist(['globalcmt_' opt.CATALOG '.mat'],'file'))
                        % catalog is broken
                        disp('Local GlobalCMT catalog is corrupted!');
                        reply=input('Attempt to re-create? Y/N [N]: ','s');
                        if(~isempty(reply) && strncmpi(reply,'y',1))
                            globalcmt_create;
                            globalcmt_update;
                            opt.CATALOG=load(['globalcmt_' opt.CATALOG]);
                            SEIZMO.GLOBALCMT.(cattype)=opt.CATALOG;
                            if(seizmoverbose); disp('Loaded Catalog!'); end
                        else
                            error('seizmo:findcmt:badCatalog',...
                                'No Valid Local GlobalCMT catalog!');
                        end
                    else
                        % no catalog...ask to make one
                        disp('No Local GlobalCMT Catalog Found!');
                        reply=input('Attempt to create? Y/N [N]: ','s');
                        if(~isempty(reply) && strncmpi(reply,'y',1))
                            globalcmt_create;
                            globalcmt_update;
                            opt.CATALOG=load(['globalcmt_' opt.CATALOG]);
                            SEIZMO.GLOBALCMT.(cattype)=opt.CATALOG;
                            if(seizmoverbose); disp('Loaded Catalog!'); end
                        else
                            error('seizmo:findcmt:badCatalog',...
                                'No Local GlobalCMT catalog!');
                        end
                    end
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

