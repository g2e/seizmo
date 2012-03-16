function []=noise_stack(indir,outdir,pair,varargin)
%NOISE_STACK    Stack noise correlation functions for noise analysis
%
%    Usage:    noise_stack(indir,outdir,pair)
%              noise_stack(indir,outdir,pair,'opt1',val,...,'optN',val)
%
%    Description:
%     NOISE_STACK(INDIR,OUTDIR,PAIR) stacks noise correlation records
%     within the directory system in INDIR and outputs the result under a
%     system of directories in OUTDIR.  PAIR indicates the correlation type
%     in shortform:
%      'ZZ','RR','RT','TR','TT'
%     or a combination of those in a cell array (see the Examples section).
%     The output stack noise correlation functions are for the entire time
%     span and are the whole NCFs not just the symmetric component.  The
%     files are named as MASTERKNAME_-_SLAVEKNAME where ??????KNAME is a
%     NET.STA.HOLE.CMP code for unique identification.  The directory
%     layout of OUTDIR is as follows:
%
%      OUTDIR
%           |
%           +-> PAIR_SPAN_XCCMP
%                             |
%                             +-> TIME
%
%     where: PAIR  is 'ZZ', 'RR', etc
%            XCCMP is 'FULL', 'SYM', etc (see option XCCMP below)
%            SPAN  is 'ALL', 'MON', etc (see option SPAN below)
%            TIME  is the time span stacked and depends on LENGTH:
%                   SPAN   |  TIME_FORMAT
%                 ---------+--------------
%                   'all'  |  yyyy.ddd.hh.mm.ss_yyyy.ddd.hh.mm.ss
%                  '3mon'  |  yyyy.MM_yyyy.MM
%                   'mon'  |  MM
%                  'yrmo'  |  yyyy.MM
%                    'wk'  |  yyyy.ddd_yyyy.ddd
%                   'day'  |  yyyy.ddd
%                   '3hr'  |  yyyy.ddd.hh.mm.ss_yyyy.ddd.hh.mm.ss
%
%     NOISE_STACK(INDIR,OUTDIR,PAIR,'OPT1',VAL,...,'OPTN',VAL) gives access
%     to several stacking and selection parameters:
%      SPAN - Controls the stacked time span.  May be of the following:
%              '3hr','day','wk','mon','yrmo','3mon','all'
%             or a cell array with a combination of those.  'mon' stacks
%             across years for a particular month (for example stacking
%             Jan 2005, 2006, & 2007).  The default is 'full'.
%      XCCMP - Controls which NCF component is output.  May be any of:
%               'full','sym','pos','neg'
%              or a cell array combination.  The default is 'full'.
%      TIMESTART - process time sections from this time on []
%      TIMEEND - process times sections before this time []
%      LATRNG - include stations in this latitude range []
%      LONRNG - include stations in this longitude range []
%      NETWORKS - include records with these network codes []
%      STATIONS - include records with these station codes []
%      STREAMS - include records with these stream codes []
%      COMPONENTS - include records with these component codes []
%      FILENAMES - limit processing to files matching this file pattern []
%      QUIETWRITE   - quietly overwrite OUTDIR (default is false)
%
%    Notes:
%     - Stack time span must exceed the timesection time span.  You cannot
%       make day or 3hr stacks from day correlations.
%
%    Header changes: SCALE (number of records in stack), DEP*
%                    Z, A, F
%
%    Examples:
%     % The typical noise processing:
%     noise_setup('some/dir','raw')
%     noise_process('raw','xc')
%     noise_stack('xc','stacks','zz')
%
%     % Stack the horizontal components:
%     noise_stack('xc','stacks',{'rr' 'rt' 'tr' 'tt'})
%
%     % Several people like to use the symmetric component:
%     noise_stack('xc','stacks','zz','xccmp','sym')
%
%    See also:  NOISE_SETUP, NOISE_PROCESS, NOISE_OVERVIEW

%     Version History:
%        June 20, 2010 - added to seizmo, fixed bug that would replace
%                        correlograms rather than combine with reversed
%                        ones, added docs
%        Sep. 21, 2010 - commented out parallel processing lines
%        Jan.  6, 2011 - update for seizmofun/solofun name change
%        Feb. 14, 2011 - no longer require 3char CMP
%        Jan. 11, 2012 - fix changeheader calls for symmetric cmp
%        Jan. 25, 2012 - inport some daydirs_stackcorr code to noise_stack
%        Jan. 27, 2012 - doc update, formalized output, option parsing
%        Jan. 31, 2012 - aggregator code written, 1st working version
%        Feb.  1, 2012 - drop type option, addrecords override & skip
%                        checking = 2x faster
%        Mar.  8, 2012 - doc update
%        Mar. 15, 2012 - parallel verbose fixes
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 15, 2012 at 11:15 GMT

% todo:
% - stack2stack
% - clear user0/1?

% check nargin
error(nargchk(3,inf,nargin));
if(nargin>=4 && ~mod(nargin,2))
    error('seizmo:noise_stack:badInput',...
        'Unpaired option/value pair given!');
end

% parse/check options
opt=noise_stack_parameters(varargin{:});

% check directories
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:noise_stack:fileNotString',...
        'INDIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:noise_stack:dirConflict',...
        ['Input Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~ischar(outdir) || ~isvector(outdir))
    error('seizmo:noise_stack:fileNotString',...
        'OUTDIR must be a string!');
end
if(exist(outdir,'file'))
    if(~exist(outdir,'dir'))
        error('seizmo:noise_stack:dirConflict',...
            'Output Directory: %s\nIs a file!',outdir);
    end
    if(~opt.QUIETWRITE)
        fprintf('Output Directory: %s\nDirectory Exists!\n',outdir);
        reply=input('Overwrite? Y/N [N]: ','s');
        if(isempty(reply) || ~strncmpi(reply,'y',1))
            disp('Not overwriting!');
            return;
        end
        disp('Overwriting!');
    end
end

% check cc
if(ischar(pair)); pair=cellstr(pair); end
if(~iscellstr(pair) || any(cellfun('prodofsize',pair)~=2))
    error('seizmo:noise_stack:badInput',...
        'PAIR must be a string like ''ZZ'' etc!');
end
pair=unique(lower(pair(:)))'; % row vector of unique strings

% directory separator
fs=filesep;

% parallel processing setup
verbose=seizmoverbose;
%matlabpool(4); % PARALLEL

% get year directories and time-section directories
dirs=xdir([indir fs]);
dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dirs
yrdir={dirs.name};
nyr=numel(yrdir);
tsdir=cell(size(yrdir));
for i=1:nyr
    dirs=xdir([indir fs yrdir{i}]);
    dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dir
    tsdir{i}={dirs.name};
end
tsdirs=[tsdir{:}]'; % all timesections
clear dirs yrdir nyr tsdir;

% convert time ranges to numeric arrays
t=char(tsdirs);
tsbgn=str2double([cellstr(t(:,1:4)) cellstr(t(:,6:8)) ...
    cellstr(t(:,10:11)) cellstr(t(:,13:14)) cellstr(t(:,16:17))]);
tsend=str2double([cellstr(t(:,19:22)) cellstr(t(:,24:26)) ...
    cellstr(t(:,28:29)) cellstr(t(:,31:32)) cellstr(t(:,34:35))]);
clear t;

% forget about any outside of user specified time limits
good=true(numel(tsdirs),1);
if(~isempty(opt.TIMESTART))
    good=good & timediff(opt.TIMESTART,tsend)>0;
end
if(~isempty(opt.TIMEEND))
    good=good & timediff(opt.TIMEEND,tsbgn)<0;
end
tsdirs=tsdirs(good); tsbgn=tsbgn(good,:); tsend=tsend(good,:);
clear good;

% get time range of timesections
[b,ib]=min(gregorian2serial(tsbgn));
[e,ie]=max(gregorian2serial(tsend));

% drop a second off of e for sanity
e=e-1/86400;

% check timesection lengths
% - adding UTC to timediff here would break the detections below
if(numel(unique(timediff(tsbgn,tsend)))>1)
    error('seizmo:noise_stack:variableTimeSectionWidth',...
        'Timesection time spans vary in INDIR!');
end
tslen=timediff(tsbgn(1,:),tsend(1,:));
if(any(ismember(opt.SPAN,'3hr')) && tslen>=10800)
    error('seizmo:noise_stack:badInput',...
        'SPAN=''3hr'' requires timesections to be less than 3 hours!');
elseif(any(ismember(opt.SPAN,'day')) && tslen>=86400)
    error('seizmo:noise_stack:badInput',...
        'SPAN=''day'' requires timesections to be less than 1 day!');
end

% detail message
if(verbose); disp('STACKING CORRELOGRAMS'); end

% turn off checking
oldseizmocheckstate=seizmocheck_state(false);
oldcheckheaderstate=checkheader_state(false);

% attempt stacking
try
    % loop over width types
    for span=opt.SPAN
        % determine time spans
        switch char(span)
            case '3hr'
                % 3hrs = 1/8th of a day
                spanbgn=serial2gregorian((fix(b*8):fix(e*8))/8,'doytime');
                spanend=serial2gregorian(((fix(b*8):fix(e*8))+1)/8,...
                    'doytime');
                spanstr=strcat(num2str(spanbgn(:,1)),'.',...
                    num2str(spanbgn(:,2),'%03d'),'.',...
                    num2str(spanbgn(:,3),'%02d'),'.',...
                    num2str(spanbgn(:,4),'%02d'),'.',...
                    num2str(spanbgn(:,5),'%02d'),'_',...
                    num2str(spanend(:,1)),'.',...
                    num2str(spanend(:,2),'%03d'),'.',...
                    num2str(spanend(:,3),'%02d'),'.',...
                    num2str(spanend(:,4),'%02d'),'.',...
                    num2str(spanend(:,5),'%02d'));
            case 'day'
                spanbgn=serial2gregorian(fix(b):fix(e),'doydate');
                spanend=serial2gregorian((fix(b):fix(e))+1,'doydate');
                spanstr=strcat(num2str(spanbgn(:,1)),'.',...
                    num2str(spanbgn(:,2),'%03d'));
            case 'wk'
                % any sun-sat with data
                % Serial Date of 2 is Sunday
                spanbgn=serial2gregorian(...
                    floor((b-2)/7)*7+2:7:floor((e-2)/7)*7+2,'doydate');
                spanend=serial2gregorian(...
                    floor((b-2)/7)*7+9:7:floor((e-2)/7)*7+9,'doydate');
                spanend1=serial2gregorian(...
                    floor((b-2)/7)*7+8:7:floor((e-2)/7)*7+8,'doydate');
                spanstr=strcat(num2str(spanbgn(:,1)),'.',...
                    num2str(spanbgn(:,2),'%03d'),'_',...
                    num2str(spanend1(:,1)),'.',...
                    num2str(spanend1(:,2),'%03d'));
                clear spanend1;
            case 'yrmo'
                btmp=serial2gregorian(b,'caldate'); btmp=btmp(1:2);
                etmp=serial2gregorian(e,'caldate'); etmp=etmp(1:2);
                spanbgn=zeros(0,3);
                spanend=spanbgn;
                for yr=btmp(1):etmp(1)
                    if(yr==btmp(1))
                        mons=btmp(2):12;
                    elseif(yr==etmp(1))
                        mons=1:etmp(2);
                    else
                        mons=1:12;
                    end
                    for mon=mons
                        spanbgn(end+1,:)=[yr mon 1];
                        spanend(end+1,:)=fixdates([yr mon+1 1]);
                    end
                end
                spanstr=strcat(num2str(spanbgn(:,1)),'.',...
                    num2str(spanbgn(:,2),'%02d'));
            case 'mon'
                % need to handle this specially
                % - mark year as nan to indicate this
                spanbgn=[nan(12,1) lind(12,1) ones(12,1)];
                spanstr=num2str(lind(12,1),'%02d');
            case '3mon'
                % only 3 month spans in which all 3 months likely have data
                btmp=serial2gregorian(b,'caldate'); btmp=btmp(1:2);
                etmp=serial2gregorian(e,'caldate'); etmp=etmp(1:2);
                spanbgn=zeros(0,3);
                spanend=spanbgn;
                spanend1=spanbgn;
                for yr=btmp(1):etmp(1)
                    if(yr==btmp(1))
                        mons=btmp(2):12;
                    elseif(yr==etmp(1))
                        mons=1:etmp(2);
                    else
                        mons=1:12;
                    end
                    for mon=mons
                        spanbgn(end+1,:)=[yr mon 1];
                        spanend(end+1,:)=fixdates([yr mon+4 1]);
                        spanend1(end+1,:)=fixdates([yr mon+3 1]);
                    end
                end
                spanbgn=spanbgn(1:end-2,:);
                spanend=spanend(1:end-2,:);
                spanstr=strcat(num2str(spanbgn(:,1)),'.',...
                    num2str(spanbgn(:,2),'%02d'),'_',...
                    num2str(spanend1(:,1)),'.',...
                    num2str(spanend1(:,2),'%02d'));
                clear spanend1;
            case 'all'
                spanbgn=tsbgn(ib,:);
                spanend=tsend(ie,:);
                spanstr=strcat(num2str(spanbgn(:,1)),'.',...
                    num2str(spanbgn(:,2),'%03d'),'.',...
                    num2str(spanbgn(:,3),'%02d'),'.',...
                    num2str(spanbgn(:,4),'%02d'),'.',...
                    num2str(spanbgn(:,5),'%02d'),'_',...
                    num2str(spanend(:,1)),'.',...
                    num2str(spanend(:,2),'%03d'),'.',...
                    num2str(spanend(:,3),'%02d'),'.',...
                    num2str(spanend(:,4),'%02d'),'.',...
                    num2str(spanend(:,5),'%02d'));
        end
        
        % loop over stack time spans
        %parfor s=1:size(spanbgn,1) % PARALLEL
        for s=1:size(spanbgn,1) % SERIAL
            % force quietness even in parfor (which resets globals)
            seizmoverbose(false);
            
            % create variables for easy access (needed for parallel)
            sbgn=spanbgn(s,:);
            send=spanend(s,:);
            
            % get timesections ENTIRELY in stack time span
            switch char(span)
                case 'mon'
                    tsbgn2=[doy2cal(tsbgn(:,1:2)) tsbgn(:,3:5)];
                    tsend2=[doy2cal(tsend(:,1:2)) tsend(:,3:5)];
                    tsend2=fixtimes([tsend2(:,1:5) tsend2(:,6)-1]); % -1s
                    in=find(tsbgn2(:,2)==sbgn(2) ...
                        & tsend2(:,2)==sbgn(2));
                otherwise
                    in=find(timediff(sbgn,tsbgn)>=0 ...
                        & timediff(send,tsend)<=0);
            end
            
            % skip in no timesections found
            if(isempty(in)); continue; end
            
            % detail message
            if(verbose); disp(['STACKING: ' spanstr(s,:)]); end
            
            % preallocation
            sdata=cell(numel(pair),1);  % stack data containers
            snamem=cell(numel(pair),1); % stack knames
            snames=cell(numel(pair),1);
            sscale=cell(numel(pair),1);
            
            % loop over timesections
            for ts=1:numel(in)
                % read in timesection data headers
                try
                    data=readheader(strcat(indir,fs,'*',fs,...
                        tsdirs(ts),fs,opt.FILENAMES));
                catch
                    % no data...
                    continue;
                end
                if(isempty(data)); continue; end
                
                % require records are correlations
                [kuser0,kuser1]=getheader(data,'kuser0','kuser1');
                xc=ismember(kuser0,{'MASTER' 'SLAVE'}) ...
                    & ismember(kuser1,{'MASTER' 'SLAVE'});
                if(~all(xc))
                    error('seizmo:noise_stack:badInput',...
                        'INDIR does not contain correlations!');
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
                
                % remove based on cmp code
                % - requiring 3char codes
                [kcmpnms,kcmpnmm]=getheader(data,'kcmpnm','kt3');
                kcmpnms=char(kcmpnms); kcmpnmm=char(kcmpnmm);
                if(size(kcmpnms,2)~=3 || size(kcmpnmm,2)~=3)
                    error('seizmo:noise_stack:badKCMPNM',...
                        'KCMPNM codes must be only 3 characters!');
                end
                kcmpnms=lower(kcmpnms(:,3)); kcmpnmm=lower(kcmpnmm(:,3));
                pidx1=ismember([kcmpnmm kcmpnms],pair);
                pidx2=ismember([kcmpnms kcmpnmm],pair);
                data=data(pidx1 | pidx2);
                
                % get some header info
                [knetwk,kstnm,khole,kcmpnm,...
                    kt0,kt1,kt2,kt3,scale]=getheader(data,...
                    'knetwk','kstnm','khole','kcmpnm',...
                    'kt0','kt1','kt2','kt3','scale');
                
                % get names
                knames=strcat(knetwk,'.',kstnm,'.',khole,'.',kcmpnm);
                knamem=strcat(kt0,'.',kt1,'.',kt2,'.',kt3);
                knames=lower(knames); knamem=lower(knamem);
                kcmpnms=char(kcmpnm); kcmpnms=lower(kcmpnms(:,3));
                kcmpnmm=char(kt3); kcmpnmm=lower(kcmpnmm(:,3));
                
                % read in data
                data=readdata(data);
                
                % multiply by scale
                scale(isnan(scale))=1;
                data=multiply(data,scale);
                
                % get reversed data
                rdata=reverse_correlations(data);
                
                % loop over pairing codes
                for p=1:numel(pair)
                    % skip if pair reversed matches a previous pair
                    % - yo dawg, data de-dup (just need to write reversed)
                    if(any(strcmp(pair{p}([2 1]),pair(1:p-1))))
                        continue;
                    end
                    
                    % skip if none have are this pair
                    pidx1=ismember([kcmpnmm kcmpnms],pair(p));
                    pidx2=ismember([kcmpnms kcmpnmm],pair(p));
                    if(~sum(pidx1 | pidx2)); continue; end
                    
                    % look up record knames in stacks
                    [tf,loc]=ismember(strcat(knamem,'.',knames),...
                        strcat(snamem{p},'.',snames{p}));
                    loc=loc(tf); % remove zeros
                    if(any(tf))
                        % those to be stacked on
                        sdata{p}(loc)=addrecords(sdata{p}(loc),data(tf),...
                            'ref','ignore');
                        sscale{p}(loc)=sscale{p}(loc)+scale(tf);
                    end
                    if(any(~tf & pidx1))
                        % those to be appended on
                        sdata{p}=[sdata{p}; data(~tf & pidx1)];
                        sscale{p}=[sscale{p}; scale(~tf & pidx1)];
                        snamem{p}=[snamem{p}; knamem(~tf & pidx1)];
                        snames{p}=[snames{p}; knames(~tf & pidx1)];
                    end
                    
                    % look up record knames reversed in stacks
                    % Note: b/c this is AFTER the above append we avoid
                    %       potential data duplication
                    [tf,loc]=ismember(strcat(knames,'.',knamem),...
                        strcat(snamem{p},'.',snames{p}));
                    loc=loc(tf); % remove zeros
                    if(any(tf))
                        % those to be stacked on
                        sdata{p}(loc)=addrecords(sdata{p}(loc),...
                            rdata(tf),'ref','ignore');
                        sscale{p}(loc)=sscale{p}(loc)+scale(tf);
                    end
                    if(any(~tf & pidx2))
                        % those to be appended on
                        sdata{p}=[sdata{p}; rdata(~tf & pidx2)];
                        sscale{p}=[sscale{p}; scale(~tf & pidx2)];
                        snamem{p}=[snamem{p}; knames(~tf & pidx2)];
                        snames{p}=[snames{p}; knamem(~tf & pidx2)];
                    end
                end
            end
            
            % data is now stacked for this time span!
            
            % loop over pairing codes
            for p=1:numel(pair)
                % de-dup
                if(any(strcmp(pair{p}([2 1]),pair(1:p-1))))
                    % which dataset?
                    pidx=find(strcmp(pair{p}([2 1]),pair(1:p-1)),1);
                    
                    % reverse
                    sdata{pidx}=reverse_correlations(sdata{pidx});
                    
                    % rename
                    sdata{pidx}=changename(sdata{pidx},...
                        'name',strcat(snames{pidx},'_-_',snamem{pidx}));
                else
                    % get span reference time
                    switch char(span)
                        case {'3hr' 'all'}
                            spanref=sbgn;
                        case {'day' 'wk' 'yrmo' '3mon'}
                            spanref=[sbgn 0 0 0];
                        case 'mon'
                            % use 1st ts's year as the year for span
                            sbgn(1)=tsbgn(in(1),1);
                            send=fixdates(sbgn+[0 1 0]);
                            spanref=[sbgn 0 0 0];
                    end
                    
                    % divide by scale to get back to an average
                    sdata{p}=divide(sdata{p},sscale{p});
                    
                    % rename
                    sdata{p}=changename(sdata{p},...
                        'name',strcat(snamem{p},'_-_',snames{p}));
                    
                    % update headers
                    sdata{p}=changeheader(sdata{p},'scale',sscale{p},...
                        'a',0,'z',spanref,'iztype','ia',...
                        'f',timediff(sbgn,send,'utc'));
                    
                    % for de-dup
                    pidx=p;
                end
                
                % xccmp
                for xc=opt.XCCMP
                    path=strcat(outdir,fs,pair(p),'_',span,'_',xc,...
                        fs,spanstr(s,:));
                    switch char(xc)
                        case 'full'
                            writeseizmo(sdata{pidx},'path',path);
                        case 'sym'
                            writeseizmo(symcmp(sdata{pidx}),'path',path);
                        case 'pos'
                            writeseizmo(cut(sdata{pidx},0),'path',path);
                        case 'neg'
                            writeseizmo(cut(reverse(sdata{pidx}),0),...
                                'path',path);
                    end
                end
            end
        end
    end
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % parallel processing takedown & fix verbosity
    %matlabpool close; % PARALLEL
    seizmoverbose(verbose);
    
    % rethrow error
    error(lasterror);
end

% toggle checking back
seizmocheck_state(oldseizmocheckstate);
checkheader_state(oldcheckheaderstate);

% parallel processing takedown & fix verbosity
%matlabpool close; % PARALLEL
seizmoverbose(verbose);

end


function [opt]=noise_stack_parameters(varargin)
% parses/checks noise_stack pv pairs

% defaults
varargin=[{'span' 'all' 'xccmp' 'full' ...
    'q' false 'ts' [] 'te' [] 'lat' [] 'lon' [] ...
    'net' [] 'sta' [] 'str' [] 'cmp' [] 'file' []} varargin];

% require option/value pairs
if(mod(nargin,2))
    error('seizmo:noise_stack:badInput',...
        'Unpaired option/value pair given!');
elseif(~iscellstr(varargin(1:2:end)))
    error('seizmo:noise_stack:badInput',...
        'Options must be specified as strings!');
end

% get user input
for i=1:2:numel(varargin)
    switch lower(varargin{i})
        case {'span' 'sp' 'spn'}
            if(isempty(varargin{i+1})); continue; end
            opt.SPAN=varargin{i+1};
        case {'xccmp' 'xc' 'xcc'}
            if(isempty(varargin{i+1})); continue; end
            opt.XCCMP=varargin{i+1};
        case {'ts' 'tstart' 'timestart'}
            opt.TIMESTART=varargin{i+1};
        case {'te' 'tend' 'timeend'}
            opt.TIMEEND=varargin{i+1};
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
        case {'q' 'qw' 'quiet' 'qwrite' 'quietwrite'}
            if(isempty(varargin{i+1})); continue; end
            opt.QUIETWRITE=varargin{i+1};
        otherwise
            error('seizmo:noise_stack:badInput',...
                'Unknown Option: %s !',varargin{i});
    end
end

% fix string options to be cellstr vectors
if(ischar(opt.SPAN)); opt.SPAN=cellstr(opt.SPAN); end
if(ischar(opt.XCCMP)); opt.XCCMP=cellstr(opt.XCCMP); end
if(ischar(opt.NETWORKS)); opt.NETWORKS=cellstr(opt.NETWORKS); end
if(ischar(opt.STATIONS)); opt.STATIONS=cellstr(opt.STATIONS); end
if(ischar(opt.STREAMS)); opt.STREAMS=cellstr(opt.STREAMS); end
if(ischar(opt.COMPONENTS)); opt.COMPONENTS=cellstr(opt.COMPONENTS); end
if(ischar(opt.FILENAMES)); opt.FILENAMES=cellstr(opt.FILENAMES); end
if(iscellstr(opt.SPAN)); opt.SPAN=unique(lower(opt.SPAN(:)))'; end
if(iscellstr(opt.XCCMP)); opt.XCCMP=unique(lower(opt.XCCMP(:)))'; end
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

% valid values
valid.SPAN={'3hr' 'day' 'wk' 'yrmo' 'mon' '3mon' 'all'};
valid.XCCMP={'full' 'pos' 'neg' 'sym'};

% check options
szs=size(opt.TIMESTART);
sze=size(opt.TIMEEND);
if(~iscellstr(opt.SPAN) || any(~ismember(opt.SPAN,valid.SPAN)))
    error('seizmo:noise_stack:badInput',...
        'SPAN option unrecognised! Use ''3hr'', ''mon'', etc!');
elseif(~iscellstr(opt.XCCMP) || any(~ismember(opt.XCCMP,valid.XCCMP)))
    error('seizmo:noise_stack:badInput',...
        'XCCMP option unrecognised! Use ''sym'', ''pos'', etc!');
elseif(~isscalar(opt.QUIETWRITE) || ~islogical(opt.QUIETWRITE))
    error('seizmo:noise_stack:badInput',...
        'QUIETWRITE flag must be a scalar logical!');
elseif(~isempty(opt.TIMESTART) && (numel(szs)>2 || szs(1)~=1 ...
        || all(szs(2)~=[2 3 5 6]) || ~isnumeric(opt.TIMESTART) ...
        || ~isreal(opt.TIMESTART)))
    error('seizmo:noise_stack:badInput',...
        'TIMESTART must be a recognized date-time vector!');
elseif(~isempty(opt.TIMEEND) && (numel(sze)>2 || sze(1)~=1 ...
        || all(sze(2)~=[2 3 5 6]) || ~isnumeric(opt.TIMEEND) ...
        || ~isreal(opt.TIMEEND)))
    error('seizmo:noise_stack:badInput',...
        'TIMEEND must be a recognized date-time vector!');
elseif(~isempty(opt.LATRNG) && (~isnumeric(opt.LATRNG) ...
        || ~isreal(opt.LATRNG) || numel(opt.LATRNG)~=2 ...
        || size(opt.LATRNG,2)~=2 || numel(size(opt.LATRNG))~=2))
    error('seizmo:noise_stack:badInput',...
        'LATRNG must be a 2 element numeric vector as [LOW HIGH]!');
elseif(~isempty(opt.LONRNG) && (~isnumeric(opt.LONRNG) ...
        || ~isreal(opt.LONRNG) || numel(opt.LONRNG)~=2 ...
        || size(opt.LONRNG,2)~=2 || numel(size(opt.LONRNG))~=2))
    error('seizmo:noise_stack:badInput',...
        'LONRNG must be a 2 element numeric vector as [LOW HIGH]!');
elseif(~isempty(opt.NETWORKS) && (~iscellstr(opt.NETWORKS)))
    error('seizmo:noise_stack:badInput',...
        'NETWORKS must be a string list of allowed network codes!');
elseif(~isempty(opt.STATIONS) && (~iscellstr(opt.STATIONS)))
    error('seizmo:noise_stack:badInput',...
        'STATIONS must be a string list of allowed station codes!');
elseif(~isempty(opt.STREAMS) && (~iscellstr(opt.STREAMS)))
    error('seizmo:noise_stack:badInput',...
        'STREAMS must be a string list of allowed stream codes!');
elseif(~isempty(opt.COMPONENTS) && (~iscellstr(opt.COMPONENTS)))
    error('seizmo:noise_stack:badInput',...
        'COMPONENTS must be a string list of allowed component codes!');
elseif(~isempty(opt.FILENAMES) && (~iscellstr(opt.FILENAMES)))
    error('seizmo:noise_stack:badInput',...
        'FILENAMES must be a string list of allowed files!');
elseif(numel(opt.FILENAMES)>1)
    error('seizmo:noise_stack:badInput',...
        'FILENAMES must be a single pattern in NOISE_STACK!');
end

% look out for xccmp options in component
if(~isempty(opt.COMPONENTS) && any(ismember(opt.COMPONENTS,...
        {'full' 'sym' 'pos' 'neg'})))
    if(~opt.QUIETWRITE)
        fprintf(...
            'Option COMPONENTS looks to have been passed XCCMP values!');
        reply=input('Problem? Y/N [N]: ','s');
        if(isempty(reply) || ~strncmpi(reply,'y',1))
            disp('Trudging ahead!');
        else
            error('seizmo:noise_stack:badInput',...
                'COMPONENTS/XCCMP option mixup!');
        end
    end
end

end


function [data]=symcmp(data)
data=solofun(changeheader(data,'b',0),@(x)[x(ceil(end/2),:); ...
    x(floor(end/2):-1:1,:)+x(ceil(end/2)+1:end,:)]);
end


function [d1]=addrecords(d1,d2,varargin)
% simple hack for speed
try
    for i=1:numel(d1); d1(i).dep=d1(i).dep+d2(i).dep; end
catch
    error('seizmo:noise_stack:badNCFs',...
        'NCFs differ in number of points! Cannot stack!');
end
end
