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
%     directories containing the records to be stacked.  This stacking is
%     "arbitrary" because it allows stacking any set of directories of
%     correlations rather than enforcing a timing constraint like in
%     NOISE_STACK.
%
%     NCFS=NOISE_STACK_ARBITRARY(INDIRS,'OPT1',VAL,...,'OPTN',VAL) gives
%     access to several selection parameters:
%      XCREVERSE  - create time-reversed correlations (false)
%      ZTRANSFORM - use Fisher's transform for stacking (true)
%      LATRNG     - include stations in this latitude range []
%      LONRNG     - include stations in this longitude range []
%      NETWORKS   - include records with these network codes []
%      STATIONS   - include records with these station codes []
%      STREAMS    - include records with these stream codes []
%      COMPONENTS - include records with these component codes []
%      FILENAMES  - limit processing to files matching this file pattern []
%
%    Notes:
%     - Header fields A, F, Z, & IZTYPE are no longer valid for the output
%       records.
%
%    Header changes: SCALE (number of records in stack), DEP*
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
%        Jan. 26, 2014 - abs path exist fix
%        May  28, 2014 - bad bugfix: averaging works correctly (yikes!),
%                        call NO_REDUNDANT_CORRELATIONS to avoid double
%                        adding due to reversed and unreversed correlations
%                        being present in a directory, improved autodetect
%                        of sac/mat i/o, added ZTRANSFORM option (defaults
%                        to true), added XCREVERSE option (defaults to
%                        false), added FILENAMES option
%        May  29, 2014 - bugfix: headers/filenames now updated, bugfix:
%                        now catches error for global variable resetting
%        July 11, 2014 - fd i/o
%        July 17, 2014 - bugfix: solofun needs func handles not strings,
%                        bugfix: convert iftype to string
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 17, 2014 at 11:15 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% directory separator
fs=filesep;

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
    end
    if(~isabspath(indirs{i})); indirs{i}=[pwd fs indirs{i}]; end
    if(~exist(indirs{i},'dir'))
        error('seizmo:noise_stack_arbitrary:dirConflict',...
            ['Input Directory: %s\n' ...
            'Does not exist (or is not a directory)!'],indirs{i});
    end
end

% parse/check options
opt=noise_stack_arbitrary_parameters(varargin{:});

% get/set verbosity
verbose=seizmoverbose(false);

% turn off checking
oldseizmocheckstate=seizmocheck_state(false);
oldcheckheaderstate=checkheader_state(false);

% for filetype checking
common_iftype=[];

% attempt stacking
try
    % detail message
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
            if(~isempty(opt.FILENAMES))
                warning('seizmo:noise_stack_arbitrary:unusedOption',...
                    'FILENAMES option ignored for MAT input!');
            end
            opt.MATIO_THIS_TIME=true;
        catch
            try
                data=readheader(strcat(indirs{ts},fs,opt.FILENAMES));
                opt.MATIO_THIS_TIME=false;
            catch
                % no data...
                continue;
            end
        end
        if(isempty(data)); continue; end
        
        % require records are correlations & all the same filetype
        [kuser0,kuser1,iftype]=getheader(data,...
            'kuser0','kuser1','iftype id');
        xc=strcmp(kuser0,'MASTER') & strcmp(kuser1,'SLAVE');
        if(~all(xc))
            error('seizmo:noise_stack_arbitrary:badInput',...
                'INDIR contains non-correlations!');
        elseif(numel(unique(iftype))~=1)
            error('seizmo:noise_stack_arbitrary:badInput',...
                'INDIR contains mixed correlation filetypes!');
        end
        
        % now check filetype is consistent across directories
        iftype=iftype{1};
        if(isempty(common_iftype))
            common_iftype=iftype;
        elseif(~strcmpi(common_iftype,iftype))
            error('seizmo:noise_stack_arbitrary:badInput',...
                'Correlations must all share the same filetype!');
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
        
        % remove redundant (reversed) correlations
        data=no_redundant_correlations(data);
        
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
        if(~opt.MATIO_THIS_TIME); data=readdata(data); end
        
        % get reversed data
        rdata=reverse_correlations(data);
        
        % convert fd to cplx
        switch lower(common_iftype)
            case 'irlim'
                data=solofun(data,@(x)complex(x(:,1),x(:,2)));
                rdata=solofun(rdata,@(x)complex(x(:,1),x(:,2)));
            case 'iamph'
                data=solofun(data,@(x)x(:,1).*exp(1j*x(:,2)));
                rdata=solofun(rdata,@(x)x(:,1).*exp(1j*x(:,2)));
        end
        
        % apply Fisher's transform
        if(opt.ZTRANS)
            data=solofun(data,@fisher);
            rdata=solofun(rdata,@fisher);
        end
        
        % multiply by scale
        % - this allows weighted stacks
        scale(isnan(scale))=1;
        data=multiply(data,scale);
        rdata=multiply(rdata,scale);
        
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
        % - don't allow autoxc this time (to avoid double add)
        [tf2,loc]=ismember(strcat(knames,'.',knamem),...
            strcat(snamem,'.',snames));
        autoxc=strcmp(knamem,knames);
        tf2=tf2 & ~autoxc;
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
        tf=~tf1 & ~tf2 & ~autoxc;
        if(opt.XCREVERSE && any(tf))
            % those to be appended on
            sdata=[sdata; rdata(tf)];
            sscale=[sscale; scale(tf)];
            snamem=[snamem; knames(tf)];
            snames=[snames; knamem(tf)];
        end
        
        % detail message
        if(verbose); print_time_left(ts,n); end
    end
    
    % divide by scale to get back to an average
    % - updates dep* stats skipped by addrecords hack
    sdata=divide(sdata,sscale);
    
    % unapply Fisher's transform
    if(opt.ZTRANS); sdata=solofun(sdata,@ifisher); end
    
    % convert cplx to fd
    % - updates dep* to not be complex
    switch lower(common_iftype)
        case 'irlim'
            sdata=solofun(sdata,@(x)[real(x),imag(x)]);
        case 'iamph'
            sdata=solofun(sdata,@(x)[abs(x),angle(x)]);
    end
    
    % rename
    sdata=changename(sdata,'name',strcat(snamem,'_-_',snames));
    
    % update headers
    sdata=changeheader(sdata,'scale',sscale);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % reset verbosity
    seizmoverbose(verbose);
    
    % rethrow error
    error(lasterror);
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
varargin=[{'lat' [] 'lon' [] 'n' [] 'st' [] 'str' [] 'cmp' [] ...
    'xcr' false 'ztrans' true 'file' []} varargin];

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
        case {'xcreverse' 'xcrev' 'xcr'}
            if(isempty(varargin{i+1})); continue; end
            opt.XCREVERSE=varargin{i+1};
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
        case {'f' 'file' 'filename' 'files' 'filenames'}
            opt.FILENAMES=varargin{i+1};
        case {'z' 'ztran' 'ztrans' 'ztransform' 'fish' 'fisher'}
            if(isempty(varargin{i+1})); continue; end
            opt.ZTRANS=varargin{i+1};
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
if(ischar(opt.FILENAMES)); opt.FILENAMES=cellstr(opt.FILENAMES); end
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
if(iscellstr(opt.FILENAMES)); opt.FILENAMES=unique(opt.FILENAMES(:)); end

% check options
if(~isscalar(opt.XCREVERSE) || ~islogical(opt.XCREVERSE))
    error('seizmo:noise_stack_arbitrary:badInput',...
        'XCREVERSE must be TRUE or FALSE!');
elseif(~isempty(opt.LATRNG) && (~isnumeric(opt.LATRNG) ...
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
elseif(~isempty(opt.FILENAMES) && (~iscellstr(opt.FILENAMES)))
    error('seizmo:noise_stack_arbitrary:badInput',...
        'FILENAMES must be a string list of allowed files!');
elseif(~isscalar(opt.ZTRANS) || ~islogical(opt.ZTRANS))
    error('seizmo:noise_stack_arbitrary:badInput',...
        'ZTRANSFORM must be TRUE or FALSE!');
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

