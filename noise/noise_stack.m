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
%            XCCMP is 'TWOWAY', 'SYMMETRIC', etc (see option XCCMP below)
%            SPAN  is 'ALL', 'MON', etc (see option SPAN below)
%            TIME  is the time span stacked and depends on SPAN:
%                   SPAN   |  TIME_FORMAT
%                 ---------+--------------
%                   'all'  |  yyyy.ddd.hh.mm.ss_yyyy.ddd.hh.mm.ss
%               '2season'  |  summer, winter
%               '4season'  |  spring, summer, autumn, winter
%                  '3mon'  |  yyyy.MM_yyyy.MM
%                   'mon'  |  MM
%                  'yrmo'  |  yyyy.MM
%                    'wk'  |  yyyy.ddd_yyyy.ddd
%                 'wkday'  |  1-SUN, 2-MON, ...
%                   'day'  |  yyyy.ddd
%                   '3hr'  |  yyyy.ddd.hh.mm.ss_yyyy.ddd.hh.mm.ss
%                   '1hr'  |  yyyy.ddd.hh.mm.ss_yyyy.ddd.hh.mm.ss
%                'tod3hr'  |  hh00 (eg. 0000, 0300, etc)
%                'tod1hr'  |  hh00 (eg. 0000, 0100, etc)
%
%     NOISE_STACK(INDIR,OUTDIR,PAIR,'OPT1',VAL,...,'OPTN',VAL) gives access
%     to several stacking and selection parameters:
%
%      SPAN - Controls the stacked time span.  May be any of the following
%             or a combination (in a cellstr array):
%              tod1hr - stack by hour of the day (across all days)
%              tod3hr - stack 3 hour blocks across all days (no overlap)
%                 1hr - stack each hour separately
%                 3hr - stack all 3 hour blocks separately
%                 day - stack each day separately
%               wkday - stack each weekday (across all weeks)
%                  wk - stack each week separately
%                 mon - stack each month (across all years)
%                yrmo - stack each month separately
%                3mon - stack each 3 month block possible (steps @ 1 month)
%             4season - stack 3 month blocks across years
%                       1-3=winter, 4-6=spring, 7-9=summer, 10-12=autumn
%             2season - stack 6 month blocks across years
%                       4-9=summer, 10-12,1-3=winter
%                 all - stack everything (DEFAULT)
%
%      XCCMP      - Controls which NCF component is output.  May be any of:
%                    'twoway','symmetric','causal','acausal'
%                   or a cell array combination.  The default is 'twoway'.
%                   XCCMP IS IGNORED FOR FREQUENCY-DOMAIN I/O.
%
%      XCREVERSE  - create time-reversed correlations (false)
%      ZTRANSFORM - use Fisher's transform for stacking (true)
%      TIMESTART  - process time sections from this time on []
%      TIMEEND    - process times sections before this time []
%      LATRNG     - include stations in this latitude range []
%      LONRNG     - include stations in this longitude range []
%      NETWORKS   - include records with these network codes []
%      STATIONS   - include records with these station codes []
%      STREAMS    - include records with these stream codes []
%      COMPONENTS - include records with these component codes []
%      FILENAMES  - limit processing to files matching this file pattern []
%      QUIETWRITE - quietly overwrite OUTDIR (default is false)
%
%    Notes:
%     - Stack time span must exceed the timesection time span.  You cannot
%       make day or 3hr stacks from day correlations.  NOISE_STACK will
%       generate an error if you attempt to do this.
%     - While data i/o is .mat files by default, SAC files are okay too.
%       Note that .mat files can be read into Matlab using the LOAD
%       function and written out as SAC files using the function
%       WRITESEIZMO.  You may also convert the inputs & outputs using
%       NOISE_MAT2SAC & NOISE_SAC2MAT.  Be aware that SAC files are ignored
%       when there is a .mat file present.  Output format is based on the
%       format of the first data file read in.
%
%    Header changes: SCALE (number of records in stack), DEP*
%                    Z, A, F
%
%    Examples:
%     % The typical noise processing:
%     noise_setup('some/dir','raw')
%     noise_process('raw','ncfs')
%     noise_stack('ncfs','stacks','zz')
%
%     % Stack the horizontal components:
%     noise_stack('ncfs','stacks',{'rr' 'rt' 'tr' 'tt'})
%
%     % Several people like to use the symmetric component:
%     noise_stack('ncfs','stacks','zz','xccmp','sym')
%
%    See also: NOISE_SETUP, NOISE_PROCESS, NOISE_OVERVIEW, STACK2STACK,
%              NOISE_STACK_ARBITRARY, NOISE_STACK_DELAZ

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
%        June 14, 2012 - xccmp value changes, make time options names more
%                        flexible
%        Aug. 22, 2012 - added 1hr, tod1hr, tod3hr, wkday, 4season, 2season
%                        span options
%        Aug. 24, 2012 - fix timespan indexing bug (span~=all were wrong)
%        Mar. 28, 2013 - matio fixes
%        Mar. 29, 2013 - autoxc bugfixes, xcreverse option
%        Apr.  3, 2013 - minor fixes and notes
%        Apr.  8, 2013 - xcreverse option is false by default (for space)
%        Sep. 24, 2013 - minor doc update
%        Jan. 26, 2014 - abs path exist fix
%        May  27, 2014 - remove MATIO option: now autodetects sac/mat i/o,
%                        added ZTRANSFORM option (defaults to true)
%        May  28, 2014 - call NO_REDUNDANT_CORRELATIONS to avoid double
%                        adding due to reversed and unreversed correlations
%                        being present in a directory, bugfix: allow
%                        multiple filenames for FILENAMES option
%        June  4, 2014 - also set t0 & t1 header fields
%        June 25, 2014 - edits for irlim fd i/o
%        July  8, 2014 - set t3-4 header fields
%        July 10, 2014 - fd is converted to complex so Fisher transform
%                        works properly (FISHER was updated), iamph ok
%        July 16, 2014 - bugfix: solofun needs func handles not strings,
%                        bugfix: convert iftype to string
%        July 21, 2014 - bugfix: use symmetric_correlations instead of
%                        local symcmp function which forgot to divide by 2
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 21, 2014 at 11:15 GMT

% todo:
% - overlap option
%   - default to 50% ?
%   - weight summation by overlap?
%   - need a time function to return overlap between ranges
%   - tod?hr, wkday, mon, ?season need special care

% check nargin
error(nargchk(3,inf,nargin));
if(nargin>=4 && ~mod(nargin,2))
    error('seizmo:noise_stack:badInput',...
        'Unpaired option/value pair given!');
end

% directory separator
fs=filesep;

% parse/check options
opt=noise_stack_parameters(varargin{:});

% check directories
if(~isstring(indir))
    error('seizmo:noise_stack:fileNotString',...
        'INDIR must be a string!');
end
if(~isabspath(indir)); indir=[pwd fs indir]; end
if(~exist(indir,'dir'))
    error('seizmo:noise_stack:dirConflict',...
        ['Input Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~isstring(outdir))
    error('seizmo:noise_stack:fileNotString',...
        'OUTDIR must be a string!');
end
if(~isabspath(outdir)); outdir=[pwd fs outdir]; end
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
% - allow those overlapping boundary
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
smaxlen={'tod1hr' 'tod3hr' '1hr' '3hr' 'day' 'wkday' ...
    'wk' 'yrmo' 'mon' '3mon' '4season' '2season' 'all'};
nmaxlen=[3601 10801 3600 10800 86400 86401 ...
    86400*7 86400*28 86400*28 86400*89 86400*90 86400*180 inf];
if(numel(unique(timediff(tsbgn,tsend)))>1)
    error('seizmo:noise_stack:variableTimeSectionWidth',...
        'Timesection time spans vary in INDIR!');
end
tslen=timediff(tsbgn(1,:),tsend(1,:));
for i=1:numel(opt.SPAN)
    if(tslen>=nmaxlen(strcmp(opt.SPAN{i},smaxlen)))
        error('seizmo:noise_stack:badInput',...
            'SPAN=''%s'' requires timesections to be less than %ds!',...
            opt.SPAN{i},nmaxlen(strcmp(opt.SPAN{i},smaxlen)));
    end
end

% detail message
if(verbose); disp('STACKING CORRELOGRAMS'); end

% turn off checking
oldseizmocheckstate=seizmocheck_state(false);
oldcheckheaderstate=checkheader_state(false);

% for autodetecting if output is .mat or SAC files
opt.MATIO=[];

% for filetype checking
common_iftype=[];

% attempt stacking
try
    % loop over width types
    for span=opt.SPAN
        % determine time spans
        switch char(span)
            case 'tod1hr'
                % need to handle this specially b/c of multiple time ranges
                spanbgn=[nan(24,2) (0:23)' zeros(24,2)];
                spanend=[nan(24,2) (1:24)' zeros(24,2)];
                spanstr=num2str((0:23)','%02d00');
            case 'tod3hr'
                % need to handle this specially b/c of multiple time ranges
                spanbgn=[nan(8,2) (0:3:21)' zeros(8,2)];
                spanend=[nan(8,2) (3:3:24)' zeros(8,2)];
                spanstr=num2str((0:3:21)','%02d00');
            case '1hr'
                % 1hr = 1/24th of a day
                spanbgn=serial2gregorian((fix(b*24):fix(e*24))/24,...
                    'doytime');
                spanend=serial2gregorian(((fix(b*24):fix(e*24))+1)/24,...
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
            case 'wkday'
                % need to handle this specially b/c of multiple time ranges
                spanbgn=[nan(7,1) (1:7)'];
                spanend=[nan(7,1) (1:7)'];
                spanstr=['1-SUN';'2-MON';'3-TUE';'4-WED';...
                    '5-THU';'6-FRI';'7-SAT'];
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
                % need to handle this specially b/c of multiple time ranges
                spanbgn=[nan(12,1) (1:12)' ones(12,1)];
                spanend=[nan(12,1) (1:12)' ones(12,1)];
                spanstr=num2str((1:12)','%02d');
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
            case '4season'
                % need to handle this specially b/c of multiple time ranges
                spanbgn=[nan(4,1) [1;4;7;10] ones(4,1)];
                spanend=[nan(4,1) [4;7;10;13] ones(4,1)];
                spanstr=['spring';'summer';'autumn';'winter'];
            case '2season'
                % need to handle this specially b/c of multiple time ranges
                spanbgn=[nan(2,1) [4;10] ones(2,1)];
                spanend=[nan(2,1) [10;4] ones(2,1)];
                spanstr=['summer';'winter'];
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
                case '2season'
                    tsbgn2=[doy2cal(tsbgn(:,1:2)) tsbgn(:,3:5)];
                    tsend2=[doy2cal(tsend(:,1:2)) tsend(:,3:5)];
                    tsend2=fixtimes([tsend2(:,1:5) tsend2(:,6)-1]); % -1s
                    % special algorithm to handle 10-12,1-3 & 4-9
                    in=find(...
                        (tsbgn2(:,2)>=sbgn(2) & tsend2(:,2)<sbgn(2)+3) ...
                        | (tsbgn2(:,2)>=send(2)-3 & tsend2(:,2)<send(2)));
                case '4season'
                    tsbgn2=[doy2cal(tsbgn(:,1:2)) tsbgn(:,3:5)];
                    tsend2=[doy2cal(tsend(:,1:2)) tsend(:,3:5)];
                    tsend2=fixtimes([tsend2(:,1:5) tsend2(:,6)-1]); % -1s
                    in=find(tsbgn2(:,2)>=sbgn(2) & tsend2(:,2)<send(2));
                case 'mon'
                    tsbgn2=[doy2cal(tsbgn(:,1:2)) tsbgn(:,3:5)];
                    tsend2=[doy2cal(tsend(:,1:2)) tsend(:,3:5)];
                    tsend2=fixtimes([tsend2(:,1:5) tsend2(:,6)-1]); % -1s
                    in=find(tsbgn2(:,2)==sbgn(2) ...
                        & tsend2(:,2)==sbgn(2));
                case {'tod1hr' 'tod3hr'}
                    tsend2=fixtimes([tsend(:,1:4) tsend(:,5)-1]); % -1s
                    in=find(tsbgn(:,3)>=sbgn(3) & tsend2(:,3)<send(3));
                case 'wkday'
                    tsend2=fixtimes([tsend(:,1:4) tsend(:,5)-1]); % -1s
                    in=find(weekday(...
                        gregorian2serial(tsbgn(:,[1 2])))==sbgn(2) ...
                        & weekday(...
                        gregorian2serial(tsend2(:,[1 2])))==sbgn(2));
                otherwise
                    in=find(timediff(sbgn,tsbgn)>=0 ...
                        & timediff(send,tsend)<=0);
            end
            
            % skip if no timesections found
            if(isempty(in)); continue; end
            
            % detail message
            nin=numel(in);
            if(verbose)
                disp(['STACKING: ' spanstr(s,:)]);
                fprintf('==> %d TIMESECTIONS IN THIS STACK\n',nin)
                print_time_left(0,nin);
            end
            
            % preallocation
            sdata=cell(numel(pair),1);  % stack data containers
            snamem=cell(numel(pair),1); % stack knames
            snames=cell(numel(pair),1);
            sscale=cell(numel(pair),1);
            
            % loop over timesections
            for ts=1:nin
                % read in timesection data headers
                try
                    data=load(strcat(indir,fs,tsdirs{in(ts)}(1:4),...
                        fs,tsdirs{in(ts)},fs,'noise_records.mat'));
                    data=data.noise_records;
                    if(~isempty(opt.FILENAMES))
                        warning('seizmo:noise_stack:unusedOption',...
                            'FILENAMES option ignored for MAT input!');
                    end
                    opt.MATIO_THIS_TIME=true;
                    if(isempty(opt.MATIO)); opt.MATIO=true; end
                catch
                    try
                        data=readheader(strcat(indir,fs,...
                            tsdirs{in(ts)}(1:4),fs,tsdirs{in(ts)},fs,...
                            opt.FILENAMES));
                        opt.MATIO_THIS_TIME=false;
                        if(isempty(opt.MATIO)); opt.MATIO=false; end
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
                    error('seizmo:noise_stack:badInput',...
                        'INDIR contains non-correlations!');
                elseif(numel(unique(iftype))~=1)
                    error('seizmo:noise_stack:badInput',...
                        'INDIR contains mixed correlation filetypes!');
                end
                
                % now check filetype is consistent across directories
                iftype=iftype{1};
                if(isempty(common_iftype))
                    common_iftype=iftype;
                elseif(~strcmpi(common_iftype,iftype))
                    error('seizmo:noise_stack:badInput',...
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
                kcmpnms=char(kcmpnm); kcmpnms=lower(kcmpnms(:,3));
                kcmpnmm=char(kt3); kcmpnmm=lower(kcmpnmm(:,3));
                
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
                % - this is also useful for stack2stack
                scale(isnan(scale))=1;
                data=multiply(data,scale);
                rdata=multiply(rdata,scale);
                
                % for debugging
                for i=1:numel(data)
                    data(i).misc.stacknames={[data(i).path data(i).name]};
                    rdata(i).misc.stacknames=...
                        {[rdata(i).path rdata(i).name]};
                end
                
                % loop over pairing codes
                for p=1:numel(pair)
                    % skip if pair reversed matches a previous pair
                    % - yo dawg, data de-dup (just need to write reversed)
                    if(any(strcmp(pair{p}([2 1]),pair(1:p-1))))
                        continue;
                    end
                    
                    % skip if none have this pair
                    pidx1=ismember([kcmpnmm kcmpnms],pair(p));
                    pidx2=ismember([kcmpnms kcmpnmm],pair(p));
                    if(~sum(pidx1 | pidx2)); continue; end
                    
                    % look up record knames in stacks
                    [tf1,loc]=ismember(strcat(knamem,'.',knames),...
                        strcat(snamem{p},'.',snames{p}));
                    loc=loc(tf1); % remove zeros
                    if(any(tf1))
                        % those to be stacked on
                        sdata{p}(loc)=addrecords(sdata{p}(loc),data(tf1),...
                            'ref','ignore');
                        sscale{p}(loc)=sscale{p}(loc)+scale(tf1);
                    end
                    
                    % look up record knames reversed in stacks
                    % - don't allow autoxc this time (to avoid double add)
                    [tf2,loc]=ismember(strcat(knames,'.',knamem),...
                        strcat(snamem{p},'.',snames{p}));
                    autoxc=strcmp(knamem,knames);
                    tf2=tf2 & ~autoxc;
                    loc=loc(tf2); % remove zeros
                    if(any(tf2))
                        % those to be stacked on
                        sdata{p}(loc)=addrecords(sdata{p}(loc),...
                            rdata(tf2),'ref','ignore');
                        sscale{p}(loc)=sscale{p}(loc)+scale(tf2);
                    end
                    
                    % append unknown records
                    tf=~tf1 & ~tf2 & pidx1;
                    if(any(tf))
                        % those to be appended on
                        sdata{p}=[sdata{p}; data(tf)];
                        sscale{p}=[sscale{p}; scale(tf)];
                        snamem{p}=[snamem{p}; knamem(tf)];
                        snames{p}=[snames{p}; knames(tf)];
                    end
                    tf=~tf1 & ~tf2 & ~autoxc & pidx2;
                    if(opt.XCREVERSE && any(tf))
                        % those to be appended on
                        sdata{p}=[sdata{p}; rdata(tf)];
                        sscale{p}=[sscale{p}; scale(tf)];
                        snamem{p}=[snamem{p}; knames(tf)];
                        snames{p}=[snames{p}; knamem(tf)];
                    end
                end
                
                % detail message
                if(verbose); print_time_left(ts,nin); end
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
                        case {'1hr' '3hr' 'all'}
                            spanref=sbgn;
                        case {'day' 'wk' 'yrmo' '3mon'}
                            spanref=[sbgn 0 0 0];
                        case {'tod1hr' 'tod3hr'}
                            % use 1st ts's info
                            sbgn(1:2)=tsbgn(in(1),1:2);
                            send(1:2)=tsbgn(in(1),1:2);
                            spanref=sbgn;
                        case 'wkday'
                            % use 1st ts's info
                            sbgn(1:2)=tsbgn(in(1),1:2);
                            send(1:2)=fixdates(sbgn+[0 1]);
                            spanref=[sbgn 0 0 0];
                        case 'mon'
                            % use 1st ts's info
                            sbgn(1)=tsbgn(in(1),1);
                            send=fixdates(sbgn+[0 1 0]);
                            spanref=[sbgn 0 0 0];
                        case '2season'
                            % use 1st ts's info
                            sbgn(1)=tsbgn(in(1),1);
                            send=fixdates(sbgn+[0 6 0]);
                            spanref=[sbgn 0 0 0];
                        case '4season'
                            % use 1st ts's info
                            sbgn(1)=tsbgn(in(1),1);
                            send=fixdates(sbgn+[0 3 0]);
                            spanref=[sbgn 0 0 0];
                    end
                    
                    % divide by scale to get back to an average
                    % - updates dep* stats skipped by addrecords hack
                    sdata{p}=divide(sdata{p},sscale{p});
                    
                    % unapply Fisher's transform
                    if(opt.ZTRANS)
                        sdata{p}=solofun(sdata{p},@ifisher);
                    end
                    
                    % convert cplx to fd
                    % - updates dep* to not be complex
                    switch lower(common_iftype)
                        case 'irlim'
                            sdata{p}=solofun(sdata{p},...
                                @(x)[real(x),imag(x)]);
                        case 'iamph'
                            sdata{p}=solofun(sdata{p},...
                                @(x)[abs(x),angle(x)]);
                    end
                    
                    % rename
                    sdata{p}=changename(sdata{p},...
                        'name',strcat(snamem{p},'_-_',snames{p}));
                    
                    % update headers
                    sdata{p}=changeheader(sdata{p},'scale',sscale{p},...
                        'a',0,'f',timediff(sbgn,send,'utc'),...
                        't0',0,'t1',timediff(sbgn,send,'utc'),...
                        't3',0,'t4',timediff(sbgn,send,'utc'),...
                        'z',spanref,'iztype','ia');
                    
                    % for de-dup
                    pidx=p;
                end
                
                % override xccmp for frequency-domain input
                xclist=opt.XCCMP;
                if(strcmpi(iftype,{'irlim' 'iamph'}))
                    xclist={'twoway'};
                end
                
                % xccmp
                for xc=xclist
                    path=char(strcat(outdir,fs,pair(p),'_',span,'_',xc,...
                        fs,spanstr(s,:)));
                    switch char(xc)
                        case 'twoway'
                            if(opt.MATIO)
                                noise_records=changepath(...
                                    sdata{pidx},'path',path); %#ok<*NASGU>
                                if(~exist(path,'dir')); mkdir(path); end
                                save([path fs 'noise_records.mat'],...
                                    'noise_records');
                                clear noise_records;
                            else
                                writeseizmo(sdata{pidx},'path',path);
                            end
                        case 'symmetric'
                            if(opt.MATIO)
                                noise_records=changepath(...
                                    symmetric_correlations(sdata{pidx}),...
                                    'path',path);
                                if(~exist(path,'dir')); mkdir(path); end
                                save([path fs 'noise_records.mat'],...
                                    'noise_records');
                                clear noise_records;
                            else
                                writeseizmo(symcmp(sdata{pidx}),...
                                    'path',path);
                            end
                        case 'causal'
                            if(opt.MATIO)
                                noise_records=changepath(...
                                    cut(sdata{pidx},0),'path',path);
                                if(~exist(path,'dir')); mkdir(path); end
                                save([path fs 'noise_records.mat'],...
                                    'noise_records');
                                clear noise_records;
                            else
                                writeseizmo(cut(sdata{pidx},0),...
                                    'path',path);
                            end
                        case 'acausal'
                            if(opt.MATIO)
                                noise_records=changepath(...
                                    cut(reverse(sdata{pidx}),0),...
                                    'path',path);
                                if(~exist(path,'dir')); mkdir(path); end
                                save([path fs 'noise_records.mat'],...
                                    'noise_records');
                                clear noise_records;
                            else
                                writeseizmo(cut(reverse(sdata{pidx}),0),...
                                    'path',path);
                            end
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
varargin=[{'span' 'all' 'xccmp' 'twoway' 'o' 50 'ztrans' true ...
    'q' false 'ts' [] 'te' [] 'lat' [] 'lon' [] 'xcr' false ...
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
        case {'xcreverse' 'xcrev' 'xcr'}
            if(isempty(varargin{i+1})); continue; end
            opt.XCREVERSE=varargin{i+1};
        case {'overlap' 'over' 'ol' 'ov' 'o' 'olap'}
            % support for this option does not exist yet!
            if(isempty(varargin{i+1})); continue; end
            opt.OVERLAP=varargin{i+1};
        case {'ts' 'tstart' 'timestart' 'startt' 'stime' 'starttime'}
            opt.TIMESTART=varargin{i+1};
        case {'te' 'tend' 'timeend' 'et' 'endtime' 'etime' 'endt'}
            opt.TIMEEND=varargin{i+1};
        case {'lat' 'la' 'lar' 'latr' 'larng' 'latitude' 'latrng'}
            opt.LATRNG=varargin{i+1};
        case {'lon' 'lo' 'lor' 'lonr' 'lorng' 'longitude' 'lonrng'}
            opt.LONRNG=varargin{i+1};
        case {'knetwk' 'n' 'net' 'netwk' 'network' 'nets' 'networks'}
            opt.NETWORKS=varargin{i+1};
        case {'kstnm' 'st' 'sta' 'stn' 'stns' 'stations' 'station'}
            if(~isempty(varargin{i+1}) && isnumeric(varargin{i+1}))
                % timestart/starttime catch
                warning('seizmo:noise_stack:badInput',...
                    'TIMESTART/STATION mixup!  Assuming TIMESTART!');
                opt.TIMESTART=varargin{i+1};
            else
                opt.STATIONS=varargin{i+1};
            end
        case {'khole' 'hole' 'holes' 'str' 'strs' 'stream' 'streams'}
            opt.STREAMS=varargin{i+1};
        case {'kcmpnm' 'cmpnm' 'cmp' 'cmps' 'component' 'components'}
            opt.COMPONENTS=varargin{i+1};
        case {'f' 'file' 'filename' 'files' 'filenames'}
            opt.FILENAMES=varargin{i+1};
        case {'q' 'qw' 'quiet' 'qwrite' 'quietwrite'}
            if(isempty(varargin{i+1})); continue; end
            opt.QUIETWRITE=varargin{i+1};
        case {'z' 'ztran' 'ztrans' 'ztransform' 'fish' 'fisher'}
            if(isempty(varargin{i+1})); continue; end
            opt.ZTRANS=varargin{i+1};
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
valid.SPAN={'tod1hr' 'tod3hr' '1hr' '3hr' 'day' 'wkday' ...
    'wk' 'yrmo' 'mon' '3mon' '4season' '2season' 'all'};
valid.XCCMP={'twoway' 'causal' 'acausal' 'symmetric'};

% check options
szs=size(opt.TIMESTART);
sze=size(opt.TIMEEND);
if(~iscellstr(opt.SPAN) || any(~ismember(opt.SPAN,valid.SPAN)))
    error('seizmo:noise_stack:badInput',...
        'SPAN option unrecognised! Use ''3hr'', ''mon'', etc!');
elseif(~iscellstr(opt.XCCMP) || any(~ismember(opt.XCCMP,valid.XCCMP)))
    error('seizmo:noise_stack:badInput',...
        'XCCMP option unrecognised! Use ''twoway'', ''symmetric'', etc!');
elseif(~isscalar(opt.XCREVERSE) || ~islogical(opt.XCREVERSE))
    error('seizmo:noise_stack:badInput',...
        'XCREVERSE must be TRUE or FALSE!');
elseif(~isscalar(opt.OVERLAP) || ~isnumeric(opt.OVERLAP) ...
        || opt.OVERLAP>100 || opt.OVERLAP<0)
    error('seizmo:noise_stack:badInput',...
        'OVERLAP must be a scalar from 0-100%%!');
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
elseif(~isscalar(opt.ZTRANS) || ~islogical(opt.ZTRANS))
    error('seizmo:noise_stack:badInput',...
        'ZTRANSFORM must be TRUE or FALSE!');
end

% look out for xccmp options in component
if(~isempty(opt.COMPONENTS) && any(ismember(opt.COMPONENTS,...
        {'twoway' 'symmetric' 'causal' 'acausal'})))
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


function [d1]=addrecords(d1,d2,varargin)
% simple hack for speed (no dep* update)
try
    for i=1:numel(d1)
        d1(i).dep=d1(i).dep+d2(i).dep;
        d1(i).misc.stacknames=...
            [d1(i).misc.stacknames; d2(i).misc.stacknames];
    end
catch
    error('seizmo:noise_stack:badNCFs',...
        'NCFs differ in number of points! Cannot stack!');
end
end
