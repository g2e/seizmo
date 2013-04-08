function [sdata]=noise_stack_arbitrary(indirs,varargin)
%NOISE_STACK_ARBITRARY    Stacks NCFS given the directories
%
%    Usage:    ncfs=noise_stack_arbitrary(indirs)
%              ncfs=noise_stack_arbitrary(indirs,'opt1',val,...,'optN',val)
%
%    Description:
%     NCFS=NOISE_STACK_ARBITRARY(INDIRS) will stack the noise correlation
%     functions (NCFs) in the input directories INDIRS and return them in
%     the SEIZMO data struct NCFS.  INDIRS must be a cell array of valid
%     directories containing the records to be stacked.
%
%     NCFS=NOISE_STACK_ARBITRARY(INDIRS,'OPT1',VAL,...,'OPTN',VAL) gives
%     access to several selection parameters:
%      LATRNG     - include stations in this latitude range []
%      LONRNG     - include stations in this longitude range []
%      NETWORKS   - include records with these network codes []
%      STATIONS   - include records with these station codes []
%      STREAMS    - include records with these stream codes []
%      COMPONENTS - include records with these component codes []
%
%    Notes:
%
%    Examples:
%     % This is equivalent to SPAN='all' for noise_stack:
%     dirs=xdir('myxc/*/');
%     dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % no . or ..
%     ncfs=noise_stack_arbitrary({dirs.name});
%
%    See also: NOISE_STACK, STACK2STACK

%     Version History:
%        Apr.  3, 2013 - initial version
%        Apr.  8, 2013 - forgot to average
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  8, 2013 at 11:15 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check indirs
if(~iscellstr(indirs))
    error('seizmo:noise_stack_arbitrary:fileNotString',...
            'INDIRS must be a cellstr array!');
end
n=numel(indirs);
for i=1:n
    if(ndims(indirs{i})~=2 || size(indirs{i},1)~=1)
        error('seizmo:noise_stack_arbitrary:fileNotString',...
            'INDIRS must all be valid strings!');
    elseif(~exist(indirs{i},'dir'))
        error('seizmo:noise_stack_arbitrary:dirConflict',...
            ['Input Directory: %s\n' ...
            'Does not exist (or is not a directory)!'],indirs{i});
    end
end

% parse/check options
opt=noise_stack_arbitrary_parameters(varargin{:});

% turn off checking
oldseizmocheckstate=seizmocheck_state(false);
oldcheckheaderstate=checkheader_state(false);

% detail message
verbose=seizmoverbose(false);
if(verbose)
    fprintf('==> %d TIMESECTIONS IN THIS STACK\n',n)
    print_time_left(0,n);
end

% preallocation
sdata=[];  % stack data
snamem=[]; % stack knames
snames=[];
sscale=[]; % stack weights

% loop over timesections
fs=filesep;
for ts=1:n
    % read in timesection data headers
    try
        data=load([indirs{ts} fs 'noise_records.mat']);
        data=data.noise_records;
        WASSAC=false;
    catch
        data=readheader(indirs{ts});
        WASSAC=true;
    end
    if(isempty(data)); continue; end
    
    % require records are correlations
    [kuser0,kuser1]=getheader(data,'kuser0','kuser1');
    xc=strcmp(kuser0,'MASTER') & strcmp(kuser1,'SLAVE');
    if(~all(xc))
        error('seizmo:noise_stack_arbitrary:badInput',...
            'INDIR contains non-correlations!');
    end
    
    % limit to the stations that the user allows
    % Note: check both fields!
    if(~isempty(opt.LATRNG))
        [stla,evla]=getheader(data,'stla','evla');
        data=data(...
            stla>=min(opt.LATRNG) & stla<=max(opt.LATRNG) ...
            & evla>=min(opt.LATRNG) & evla<=max(opt.LATRNG));
        if(isempty(data)); continue; end
    end
    if(~isempty(opt.LONRNG))
        [stlo,evlo]=getheader(data,'stlo','evlo');
        data=data(...
            stlo>=min(opt.LONRNG) & stlo<=max(opt.LONRNG) ...
            & evlo>=min(opt.LONRNG) & evlo<=max(opt.LONRNG));
        if(isempty(data)); continue; end
    end
    if(~isempty(opt.NETWORKS))
        [knetwk1,knetwk2]=getheader(data,'knetwk','kt0');
        data=data(ismember(lower(knetwk1),opt.NETWORKS) ...
            & ismember(lower(knetwk2),opt.NETWORKS));
        if(isempty(data)); continue; end
    end
    if(~isempty(opt.STATIONS))
        [kstnm1,kstnm2]=getheader(data,'kstnm','kt1');
        data=data(ismember(lower(kstnm1),opt.STATIONS) ...
            & ismember(lower(kstnm2),opt.STATIONS));
        if(isempty(data)); continue; end
    end
    if(~isempty(opt.STREAMS))
        [khole1,khole2]=getheader(data,'khole','kt2');
        data=data(ismember(lower(khole1),opt.STREAMS) ...
            & ismember(lower(khole2),opt.STREAMS));
        if(isempty(data)); continue; end
    end
    if(~isempty(opt.COMPONENTS))
        [kcmpnm1,kcmpnm2]=getheader(data,'kcmpnm','kt3');
        data=data(ismember(lower(kcmpnm1),opt.COMPONENTS) ...
            & ismember(lower(kcmpnm2),opt.COMPONENTS));
        if(isempty(data)); continue; end
    end
    
    % get some header info
    [knetwk,kstnm,khole,kcmpnm,...
        kt0,kt1,kt2,kt3,scale]=getheader(data,...
        'knetwk','kstnm','khole','kcmpnm',...
        'kt0','kt1','kt2','kt3','scale');
    
    % get names
    knames=strcat(knetwk,'.',kstnm,'.',khole,'.',kcmpnm);
    knamem=strcat(kt0,'.',kt1,'.',kt2,'.',kt3);
    knames=lower(knames); knamem=lower(knamem);
    
    % read in data (if not MATFILE input)
    if(WASSAC); data=readdata(data); end
    
    % multiply by scale
    % - this allows weighted stacks
    scale(isnan(scale))=1;
    data=multiply(data,scale);
    
    % get reversed data
    rdata=reverse_correlations(data);
    
    % for debugging
    for i=1:numel(data)
        data(i).misc.stacknames={[data(i).path data(i).name]};
        rdata(i).misc.stacknames={[rdata(i).path rdata(i).name]};
    end
    
    % look up record knames in stacks
    [tf1,loc]=ismember(strcat(knamem,'.',knames),...
        strcat(snamem,'.',snames));
    loc=loc(tf1); % remove zeros
    if(any(tf1))
        % those to be stacked on
        sdata(loc)=addrecords(sdata(loc),data(tf1),'ref','ignore');
        sscale(loc)=sscale(loc)+scale(tf1);
    end
    
    % look up record knames reversed in stacks
    % - don't allow adding unreversed & reversed
    [tf2,loc]=ismember(strcat(knames,'.',knamem),...
        strcat(snamem,'.',snames));
    tf2=tf2 & ~tf1;
    loc=loc(tf2); % remove zeros
    if(any(tf2))
        % those to be stacked on
        sdata(loc)=addrecords(sdata(loc),rdata(tf2),'ref','ignore');
        sscale(loc)=sscale(loc)+scale(tf2);
    end
    
    % append unknown records
    tf=~tf1 & ~tf2;
    if(any(tf))
        % those to be appended on
        sdata=[sdata; data(tf)];
        sscale=[sscale; scale(tf)];
        snamem=[snamem; knamem(tf)];
        snames=[snames; knames(tf)];
    end
    
    % divide by scale to get back to an average
    % - updates dep* stats skipped by addrecords hack
    sdata=divide(sdata,sscale);
    
    % detail message
    if(verbose); print_time_left(ts,n); end
end

% toggle checking back
seizmocheck_state(oldseizmocheckstate);
checkheader_state(oldcheckheaderstate);

% reset verbosity
seizmoverbose(verbose);

end
function [opt]=noise_stack_arbitrary_parameters(varargin)
% parses/checks noise_stack_arbitrary pv pairs

% defaults
varargin=[{'lat' [] 'lon' [] 'n' [] 'st' [] 'str' [] 'cmp' []} varargin];

% require option/value pairs
if(mod(nargin,2))
    error('seizmo:noise_stack_arbitrary:badInput',...
        'Unpaired option/value pair given!');
elseif(~iscellstr(varargin(1:2:end)))
    error('seizmo:noise_stack_arbitrary:badInput',...
        'Options must be specified as strings!');
end

% get user input
for i=1:2:numel(varargin)
    switch lower(varargin{i})
        case {'lat' 'la' 'lar' 'latr' 'larng' 'latitude' 'latrng'}
            opt.LATRNG=varargin{i+1};
        case {'lon' 'lo' 'lor' 'lonr' 'lorng' 'longitude' 'lonrng'}
            opt.LONRNG=varargin{i+1};
        case {'knetwk' 'n' 'net' 'netwk' 'network' 'nets' 'networks'}
            opt.NETWORKS=varargin{i+1};
        case {'kstnm' 'st' 'sta' 'stn' 'stns' 'stations' 'station'}
            opt.STATIONS=varargin{i+1};
        case {'khole' 'hole' 'holes' 'str' 'strs' 'stream' 'streams'}
            opt.STREAMS=varargin{i+1};
        case {'kcmpnm' 'cmpnm' 'cmp' 'cmps' 'component' 'components'}
            opt.COMPONENTS=varargin{i+1};
        otherwise
            error('seizmo:noise_stack_arbitrary:badInput',...
                'Unknown Option: %s !',varargin{i});
    end
end

% fix string options to be cellstr vectors
if(ischar(opt.NETWORKS)); opt.NETWORKS=cellstr(opt.NETWORKS); end
if(ischar(opt.STATIONS)); opt.STATIONS=cellstr(opt.STATIONS); end
if(ischar(opt.STREAMS)); opt.STREAMS=cellstr(opt.STREAMS); end
if(ischar(opt.COMPONENTS)); opt.COMPONENTS=cellstr(opt.COMPONENTS); end
if(iscellstr(opt.NETWORKS))
    opt.NETWORKS=unique(lower(opt.NETWORKS(:)));
end
if(iscellstr(opt.STATIONS))
    opt.STATIONS=unique(lower(opt.STATIONS(:)));
end
if(iscellstr(opt.STREAMS)); opt.STREAMS=unique(lower(opt.STREAMS(:))); end
if(iscellstr(opt.COMPONENTS))
    opt.COMPONENTS=unique(lower(opt.COMPONENTS(:)));
end

% check options
if(~isempty(opt.LATRNG) && (~isnumeric(opt.LATRNG) ...
        || ~isreal(opt.LATRNG) || numel(opt.LATRNG)~=2 ...
        || size(opt.LATRNG,2)~=2 || numel(size(opt.LATRNG))~=2))
    error('seizmo:noise_stack_arbitrary:badInput',...
        'LATRNG must be a 2 element numeric vector as [LOW HIGH]!');
elseif(~isempty(opt.LONRNG) && (~isnumeric(opt.LONRNG) ...
        || ~isreal(opt.LONRNG) || numel(opt.LONRNG)~=2 ...
        || size(opt.LONRNG,2)~=2 || numel(size(opt.LONRNG))~=2))
    error('seizmo:noise_stack_arbitrary:badInput',...
        'LONRNG must be a 2 element numeric vector as [LOW HIGH]!');
elseif(~isempty(opt.NETWORKS) && (~iscellstr(opt.NETWORKS)))
    error('seizmo:noise_stack_arbitrary:badInput',...
        'NETWORKS must be a string list of allowed network codes!');
elseif(~isempty(opt.STATIONS) && (~iscellstr(opt.STATIONS)))
    error('seizmo:noise_stack_arbitrary:badInput',...
        'STATIONS must be a string list of allowed station codes!');
elseif(~isempty(opt.STREAMS) && (~iscellstr(opt.STREAMS)))
    error('seizmo:noise_stack_arbitrary:badInput',...
        'STREAMS must be a string list of allowed stream codes!');
elseif(~isempty(opt.COMPONENTS) && (~iscellstr(opt.COMPONENTS)))
    error('seizmo:noise_stack_arbitrary:badInput',...
        'COMPONENTS must be a string list of allowed component codes!');
end

end


function [d1]=addrecords(d1,d2,varargin)
% simple hack for speed (no dep* update)
try
    for i=1:numel(d1)
        d1(i).dep=d1(i).dep+d2(i).dep;
        d1(i).misc.stacknames=...
            [d1(i).misc.stacknames; d2(i).misc.stacknames];
    end
catch
    error('seizmo:noise_stack_arbitrary:badNCFs',...
        'NCFs differ in number of points! Cannot stack!');
end
end

