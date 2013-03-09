function [cmts]=findcmts(varargin)
%FINDCMTS    Returns Global CMTs in specified time & position ranges
%
%    Usage:    cmt=findcmts('option1',value1,...,'optionN',valueN)
%
%    Description:
%     CMT=FINDCMTS('OPTION1',VALUE1,...,'OPTIONN',VALUEN) finds the Global
%     CMT(s) (www.globalcmt.org) using the parameters provided.  The
%     options are:
%      TYPE          -- location type: 'centroid' or 'hypocenter' (default)
%      CATALOG       -- CMT catalog: 'Full', 'Quick', 'Both' (Default)
%      STARTTIME     -- Start of time range (see Notes for formats)
%      ENDTIME       -- End of time range (see Notes for formats)
%      NUMDAYS       -- Number of days after start time (resets ENDTIME)
%      NUMSECS       -- Number of seconds after start time (resets ENDTIME)
%      LATRANGE      -- Latitude range in degrees ([lo hi])
%      LONRANGE      -- Longitude range in degrees ([lo hi])
%      DEPRANGE      -- Depth range in kilometers ([lo hi])
%      M(W,S,B)RANGE -- Magnitude range ([lo hi])
%      HALFDURATION  -- half duration range in seconds ([lo hi])
%      CENTROIDSHIFT -- centroid time shift range in seconds ([lo hi])
%      TAXISPLUNGE   -- tension axis plunge in degrees ([lo hi])
%      NAXISPLUNGE   -- null axis plunge in degrees  ([lo hi])
%      NAME          -- limit CMTs to those with names matching the given
%                       regular expression(s) (see REGEXP for details)
%     
%     Note that NUMDAYS & NUMSECS set ENDTIME indirectly so only the last
%     call to one of the 3 will be honored.  So if you set NUMDAYS to 3 and
%     then NUMSECS to 3600, the search range is only an hour NOT 3 days +
%     1 hour.  The defaults are set to match the GlobalCMT site except for
%     the STARTTIME & ENDTIME -- they are set to return all as I find that
%     more useful than a default that returns results from the first day of
%     1976.
%
%    Notes:
%     - Longitude range may cross the dateline.
%     - Time and location constraints are placed on the centroid time and
%       location not the hypocenter.  Use the TYPE option to change that.
%     - TIME formats: [YEAR JULDAY] (start of day)
%                     [YEAR MONTH CALDAY] (start of day)
%                     [YEAR JULDAY HOUR MINUTE SECONDS]
%                     [YEAR MONTH CALDAY HOUR MINUTE SECONDS]
%     - CATALOG may be an NDK struct for searching a custom CMT catalog.
%     - GlobalCMT catalogs are cached under global SEIZMO.GLOBALCMT
%
%    Examples:
%     % Return only events with the following names:
%     cmtnames={'C031995G' 'M051695F' 'C062495A' 'C062995D'};
%     cmts=findcmts('names',cmtnames);
%
%     % Now show info for only C031995G:
%     ssidx(cmts,strcmp(cmts.name,'C031995G'))
%
%    See also: FINDCMT, SETEVENT, READNDK, SSIDX, GLOBALCMT_UPDATE

%     Version History:
%        Aug.  2, 2010 - initial version
%        Aug. 10, 2010 - add name option
%        Jan. 26, 2011 - support cellstr input into .name, altered
%                        STARTTIME & ENDTIME to [1962 1] & now, added
%                        examples
%        Feb. 11, 2011 - fix name/nullaxis option conflict
%        Mar. 11, 2011 - fix null axis bug, add lonrng/latrng, fix name bug
%        June 11, 2011 - minor defaults fix
%        Feb. 29, 2012 - doc update, type option & switch to hypo
%        Mar.  1, 2012 - auto-creation of catalogs
%        Feb. 26, 2013 - sorted by time (slower but nice)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 26, 2013 at 21:30 GMT

% todo:

% check optional inputs
opt=check_findcmts_options(varargin{:});

% get moment magnitude
mw=(2/3).*(log10(opt.CATALOG.scalarmoment.*10.^opt.CATALOG.exponent)-16.1);

% get indexing
switch lower(opt.TYPE)
    case {'centroid' 'cmt' 'c'}
        idx=(timediff(opt.STARTTIME,[opt.CATALOG.year opt.CATALOG.month ...
            opt.CATALOG.day opt.CATALOG.hour opt.CATALOG.minute ...
            opt.CATALOG.seconds+opt.CATALOG.centroidtime])>=0) ...
            & (timediff(opt.ENDTIME,[opt.CATALOG.year opt.CATALOG.month ...
            opt.CATALOG.day opt.CATALOG.hour opt.CATALOG.minute ...
            opt.CATALOG.seconds+opt.CATALOG.centroidtime])<=0) ...
            & (opt.LATRANGE(1)<=opt.CATALOG.centroidlat) ...
            & (opt.LATRANGE(2)>=opt.CATALOG.centroidlat) ...
            & inlonrng(opt.LONRANGE,opt.CATALOG.centroidlon) ...
            & (opt.DEPRANGE(1)<=opt.CATALOG.centroiddep) ...
            & (opt.DEPRANGE(2)>=opt.CATALOG.centroiddep) ...
            & (opt.MWRANGE(1)<=mw) ...
            & (opt.MWRANGE(2)>=mw) ...
            & (opt.MSRANGE(1)<=opt.CATALOG.ms) ...
            & (opt.MSRANGE(2)>=opt.CATALOG.ms) ...
            & (opt.MBRANGE(1)<=opt.CATALOG.mb) ...
            & (opt.MBRANGE(2)>=opt.CATALOG.mb) ...
            & (opt.HALFDURATION(1)<=opt.CATALOG.srcfuncdur) ...
            & (opt.HALFDURATION(2)>=opt.CATALOG.srcfuncdur) ...
            & (opt.CENTROIDSHIFT(1)<=opt.CATALOG.centroidtime) ...
            & (opt.CENTROIDSHIFT(2)>=opt.CATALOG.centroidtime) ...
            & (opt.TAXISPLUNGE(1)<=opt.CATALOG.plunge1) ...
            & (opt.TAXISPLUNGE(2)>=opt.CATALOG.plunge1) ...
            & (opt.NAXISPLUNGE(1)<=opt.CATALOG.plunge2) ...
            & (opt.NAXISPLUNGE(2)>=opt.CATALOG.plunge2);
    case {'hypocenter' 'hypo' 'h'}
        idx=(timediff(opt.STARTTIME,[opt.CATALOG.year opt.CATALOG.month ...
            opt.CATALOG.day opt.CATALOG.hour opt.CATALOG.minute ...
            opt.CATALOG.seconds])>=0) ...
            & (timediff(opt.ENDTIME,[opt.CATALOG.year opt.CATALOG.month ...
            opt.CATALOG.day opt.CATALOG.hour opt.CATALOG.minute ...
            opt.CATALOG.seconds])<=0) ...
            & (opt.LATRANGE(1)<=opt.CATALOG.latitude) ...
            & (opt.LATRANGE(2)>=opt.CATALOG.latitude) ...
            & inlonrng(opt.LONRANGE,opt.CATALOG.longitude) ...
            & (opt.DEPRANGE(1)<=opt.CATALOG.depth) ...
            & (opt.DEPRANGE(2)>=opt.CATALOG.depth) ...
            & (opt.MWRANGE(1)<=mw) ...
            & (opt.MWRANGE(2)>=mw) ...
            & (opt.MSRANGE(1)<=opt.CATALOG.ms) ...
            & (opt.MSRANGE(2)>=opt.CATALOG.ms) ...
            & (opt.MBRANGE(1)<=opt.CATALOG.mb) ...
            & (opt.MBRANGE(2)>=opt.CATALOG.mb) ...
            & (opt.HALFDURATION(1)<=opt.CATALOG.srcfuncdur) ...
            & (opt.HALFDURATION(2)>=opt.CATALOG.srcfuncdur) ...
            & (opt.CENTROIDSHIFT(1)<=opt.CATALOG.centroidtime) ...
            & (opt.CENTROIDSHIFT(2)>=opt.CATALOG.centroidtime) ...
            & (opt.TAXISPLUNGE(1)<=opt.CATALOG.plunge1) ...
            & (opt.TAXISPLUNGE(2)>=opt.CATALOG.plunge1) ...
            & (opt.NAXISPLUNGE(1)<=opt.CATALOG.plunge2) ...
            & (opt.NAXISPLUNGE(2)>=opt.CATALOG.plunge2);
end

% name idx
if(isfield(opt,'NAME') && ~isempty(opt.NAME))
    idx0=false(size(idx));
    for i=1:numel(opt.NAME)
        idx0=idx0 | ~cellfun('isempty',...
            regexp(opt.CATALOG.name,opt.NAME{i}));
    end
    idx=idx & idx0;
end

% extract cmts
cmts=ssidx(opt.CATALOG,idx);

% sort by time
[idx,idx]=sort(datenum([cmts.year cmts.month cmts.day cmts.hour ...
    cmts.minute cmts.seconds+cmts.centroidtime]));
fields=fieldnames(cmts);
for i=1:numel(fields); cmts.(fields{i})=cmts.(fields{i})(idx); end

end


function [opt]=check_findcmts_options(varargin)

% global SEIZMO access for catalog caching
global SEIZMO

% require option/value pairings
if(mod(nargin,2))
    error('seizmo:findcmts:badInput',...
        'Unpaired OPTION/VALUE!');
end

% defaults
varargin=[{'st' [1962 1] 'et' datevec(now) 'cat' 'both' 'mw' [0 10] ...
    'ms' [0 10] 'mb' [0 10] 'lat' [-90 90] 'lon' [-180 180] ...
    'dep' [0 1000] 'hd' [-9999 9999] 'cs' [-9999 9999] ...
    'tax' [0 90] 'nax' [0 90]} 'type' 'hypo' varargin];

% require options are strings
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:findcmts:badInput',...
        'OPTION must be a string!');
end

% valid strings
valid.TYPE={'hypocenter' 'hypo' 'h' 'centroid' 'cmt' 'c'};
valid.CATALOG={'FULL' 'QUICK' 'BOTH'};
valid.TIMES={'[YEAR JDAY]' '[YEAR MONTH CDAY]' ...
    '[YEAR JDAY HOUR MINUTE SECONDS]' ...
    '[YEAR MONTH CDAY HOUR MINUTE SECONDS]'};
valid.NDKFIELDS={'scalarmoment' 'exponent' 'year' 'month' 'day' 'hour' ...
    'minute' 'seconds' 'centroidtime' 'centroidlat' 'centroidlon' ...
    'centroiddep' 'srcfuncdur' 'ms' 'mb' 'plunge1' 'plunge2' ...
    'latitude' 'longitude' 'depth'};

% loop over options
for i=1:2:numel(varargin)
    switch lower(varargin{i})
        case {'names' 'name' 'n'}
            if(ischar(varargin{i+1}))
                varargin{i+1}=cellstr(varargin{i+1});
            end
            if(~isempty(varargin{i+1}) && ~iscellstr(varargin{i+1}))
                error('seizmo:findcmts:badInput',...
                    'NAME must be a string or cell string array!');
            end
            opt.NAME=varargin{i+1};
        case {'type' 'ty'}
            if(ischar(varargin{i+1}) ...
                    && ismember(lower(varargin{i+1}),valid.TYPE))
                opt.TYPE=varargin{i+1};
            else
                error('seizmo:findcmts:badInput',...
                    ['TYPE must be one of the following:\n' ...
                    sprintf('''%s'' ',valid.TYPE{:})]);
            end
        case {'starttime' 'st'}
            % check
            sz=size(varargin{i+1});
            if(~isreal(varargin{i+1}) || numel(sz)>2 ...
                    || sz(1)~=1 || ~any(sz(2)==[2 3 5 6]))
                error('seizmo:findcmts:badInput',...
                    ['STARTTIME must have one of the ' ...
                    'following datetime formats:\n' ...
                    sprintf('%s\n',valid.TIMES)]);
            end
            opt.STARTTIME=varargin{i+1};
        case {'endtime' 'et'}
            % delay interpretation until after all other options
            sz=size(varargin{i+1});
            if(~isreal(varargin{i+1}) || numel(sz)>2 ...
                    || sz(1)~=1 || ~any(sz(2)==[2 3 5 6]))
                error('seizmo:findcmts:badInput',...
                    ['ENDTIME must have one of the ' ...
                    'following datetime formats:\n' ...
                    sprintf('%s\n',valid.TIMES)]);
            end
            et=varargin(i:i+1);
        case {'numdays' 'nd'}
            % delay interpretation until after all other options
            if(~isreal(varargin{i+1}) || ~isscalar(varargin{i+1}) ...
                    || varargin{i+1}~=fix(varargin{i+1}) ...
                    || varargin{i+1}<0)
                error('seizmo:findcmts:badInput',...
                    'NUMDAYS must be a positive scalar integer!');
            end
            et=varargin(i:i+1);
        case {'numsecs' 'ns'}
            % delay interpretation until after all other options
            if(~isreal(varargin{i+1}) || ~isscalar(varargin{i+1}) ...
                    || varargin{i+1}<0)
                error('seizmo:findcmts:badInput',...
                    'NUMSECS must be a real-valued positive scalar!');
            end
            et=varargin(i:i+1);
        case {'latrange' 'latrng' 'lat' 'la'}
            % check
            if(~isreal(varargin{i+1}) ...
                    || ~isequal(size(varargin{i+1}),[1 2]) ...
                    || varargin{i+1}(1)>varargin{i+1}(2))
                error('seizmo:findcmts:badInput',...
                    'LATRANGE must be in the format [LO HI]!');
            end
            opt.LATRANGE=varargin{i+1};
        case {'lonrange' 'lonrng' 'lon' 'lo'}
            % check
            if(~isreal(varargin{i+1}) ...
                    || ~isequal(size(varargin{i+1}),[1 2]) ...
                    || varargin{i+1}(1)>varargin{i+1}(2))
                error('seizmo:findcmts:badInput',...
                    'LONRANGE must be in the format [LO HI]!');
            end
            opt.LONRANGE=varargin{i+1};
        case {'deprange' 'dep' 'd' 'z'}
            % check
            if(~isreal(varargin{i+1}) ...
                    || ~isequal(size(varargin{i+1}),[1 2]) ...
                    || varargin{i+1}(1)>varargin{i+1}(2))
                error('seizmo:findcmts:badInput',...
                    'DEPRANGE must be in the format [LO HI]!');
            end
            % auto-fix meter depths
            if(any(varargin{i+1}>1000))
                warning('seizmo:findcmts:badDepth',...
                    ['Assuming DEPTH>1000 is in meters.\n' ...
                    'Converting depths to kilometers!']);
                varargin{i+1}=varargin{i+1}/1000;
            end
            opt.DEPRANGE=varargin{i+1};
        case {'catalog' 'cat' 'c'}
            % check
            if(isstruct(varargin{i+1}) && all(ismember(valid.NDKFIELDS,...
                    fieldnames(varargin{i+1}))))
                opt.CATALOG=varargin{i+1};
            elseif(ischar(varargin{i+1}) ...
                    && ismember(upper(varargin{i+1}),valid.CATALOG))
                opt.CATALOG=varargin{i+1};
            else
                error('seizmo:findcmts:badInput',...
                    ['CATALOG must be a ndk-struct or:\n' ...
                    sprintf('''%s'' ',valid.CATALOG{:})]);
            end
        case {'mwrange' 'mw'}
            % check
            if(~isreal(varargin{i+1}) ...
                    || ~isequal(size(varargin{i+1}),[1 2]) ...
                    || varargin{i+1}(1)>varargin{i+1}(2))
                error('seizmo:findcmts:badInput',...
                    'MWRANGE must be in the format [LO HI]!');
            end
            opt.MWRANGE=varargin{i+1};
        case {'msrange' 'ms'}
            % check
            if(~isreal(varargin{i+1}) ...
                    || ~isequal(size(varargin{i+1}),[1 2]) ...
                    || varargin{i+1}(1)>varargin{i+1}(2))
                error('seizmo:findcmts:badInput',...
                    'MSRANGE must be in the format [LO HI]!');
            end
            opt.MSRANGE=varargin{i+1};
        case {'mbrange' 'mb'}
            % check
            if(~isreal(varargin{i+1}) ...
                    || ~isequal(size(varargin{i+1}),[1 2]) ...
                    || varargin{i+1}(1)>varargin{i+1}(2))
                error('seizmo:findcmts:badInput',...
                    'MBRANGE must be in the format [LO HI]!');
            end
            opt.MBRANGE=varargin{i+1};
        case {'halfduration' 'hd'}
            % check
            if(~isreal(varargin{i+1}) ...
                    || ~isequal(size(varargin{i+1}),[1 2]) ...
                    || varargin{i+1}(1)>varargin{i+1}(2))
                error('seizmo:findcmts:badInput',...
                    'HALFDURATION must be in the format [LO HI]!');
            end
            opt.HALFDURATION=varargin{i+1};
        case {'centroidshift' 'cen' 'cs'}
            % check
            if(~isreal(varargin{i+1}) ...
                    || ~isequal(size(varargin{i+1}),[1 2]) ...
                    || varargin{i+1}(1)>varargin{i+1}(2))
                error('seizmo:findcmts:badInput',...
                    'CENTROIDSHIFT must be in the format [LO HI]!');
            end
            opt.CENTROIDSHIFT=varargin{i+1};
        case {'taxisplunge' 'tension' 'tax' 'ta'}
            % check
            if(~isreal(varargin{i+1}) ...
                    || ~isequal(size(varargin{i+1}),[1 2]) ...
                    || varargin{i+1}(1)>varargin{i+1}(2))
                error('seizmo:findcmts:badInput',...
                    'TAXISPLUNGE must be in the format [LO HI]!');
            end
            opt.TAXISPLUNGE=varargin{i+1};
        case {'naxisplunge' 'null' 'nax' 'na'}
            % check
            if(~isreal(varargin{i+1}) ...
                    || ~isequal(size(varargin{i+1}),[1 2]) ...
                    || varargin{i+1}(1)>varargin{i+1}(2))
                error('seizmo:findcmts:badInput',...
                    'NAXISPLUNGE must be in the format [LO HI]!');
            end
            opt.NAXISPLUNGE=varargin{i+1};
        otherwise
            error('seizmo:findcmts:badInput',...
                'Unknown Option: %s !',varargin{i});
    end
end

% now get get endtime
switch lower(et{1})
    case {'endtime' 'et'}
        opt.ENDTIME=et{2};
    case {'numdays' 'nd'}
        switch numel(opt.STARTTIME)
            case 2
                opt.ENDTIME=fixdates(opt.STARTTIME+[0 et{2}]);
            case 3
                opt.ENDTIME=fixdates(opt.STARTTIME+[0 0 et{2}]);
            case 5
                opt.ENDTIME=fixtimes(opt.STARTTIME+[0 et{2} 0 0 0]);
            case 6
                opt.ENDTIME=fixtimes(opt.STARTTIME+[0 0 et{2} 0 0 0]);
        end
    case {'numsecs' 'ns'}
        switch numel(opt.STARTTIME)
            case 2
                opt.ENDTIME=fixtimes([opt.STARTTIME 0 0 0] ...
                    +[0 0 0 0 et{2}]);
            case 3
                opt.ENDTIME=fixtimes([opt.STARTTIME 0 0 0] ...
                    +[0 0 0 0 0 et{2}]);
            case 5
                opt.ENDTIME=fixtimes(opt.STARTTIME+[0 et{2} 0 0 0]);
            case 6
                opt.ENDTIME=fixtimes(opt.STARTTIME+[0 0 et{2} 0 0 0]);
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
                            error('seizmo:findcmts:badCatalog',...
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
                            error('seizmo:findcmts:badCatalog',...
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
                            error('seizmo:findcmts:badCatalog',...
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
                            error('seizmo:findcmts:badCatalog',...
                                'No Local GlobalCMT Quick catalog!');
                        end
                    end
                end
            end
            if(isstruct(tmp1) && ~all(ismember(valid.NDKFIELDS,...
                    fieldnames(tmp1))))
                error('seizmo:findcmts:badCatalog',...
                    'GlobalCMT catalog is corrupt!');
            end
            if(isstruct(tmp2) && ~all(ismember(valid.NDKFIELDS,...
                    fieldnames(tmp2))))
                error('seizmo:findcmts:badCatalog',...
                    'GlobalCMT catalog is corrupt!');
            end
            fields=fieldnames(tmp1);
            if(~isequal(fields,fieldnames(tmp2)))
                error('seizmo:findcmts:badCatalog',...
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
                            error('seizmo:findcmts:badCatalog',...
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
                            error('seizmo:findcmts:badCatalog',...
                                'No Local GlobalCMT catalog!');
                        end
                    end
                end
            end
            if(isstruct(opt.CATALOG) && ~all(ismember(valid.NDKFIELDS,...
                    fieldnames(opt.CATALOG))))
                error('seizmo:findcmts:badCatalog',...
                    'GlobalCMT catalog is corrupt!');
            end
    end
end

end
