function []=stack2stack(stackdir,newspan,varargin)
%STACK2STACK    Stack noise correlation function stacks
%
%    Usage:    stack2stack(stackdir,newspan)
%              stack2stack(stackdir,newspan,'opt1',val,...,'optN',val)
%
%    Description:
%     STACK2STACK(STACKDIR,NEWSPAN) stacks noise correlation stacks within
%     the directory STACKDIR and outputs the resulting stacks at the same
%     directory level with the name determined by STACKDIR & NEWSPAN.  The
%     directory layout of OUTDIR (the directory above STACKDIR and created
%     by NOISE_STACK) is as follows:
%
%      OUTDIR
%           |
%           +-> PAIR_OLDSPAN_XCCMP
%           +-> PAIR_NEWSPAN_XCCMP
%                                |
%                                +-> TIME
%
%     where TIME has a format dependent on NEWSPAN:
%                  NEWSPAN |  TIME
%                 ---------+--------------
%                   'all'  |  yyyy.ddd.hh.mm.ss_yyyy.ddd.hh.mm.ss
%               '2season'  |  summer, winter
%               '4season'  |  spring, summer, autumn, winter
%                   'mon'  |  MM
%                  '3mon'  |  yyyy.MM_yyyy.MM
%                  'yrmo'  |  yyyy.MM
%                    'wk'  |  yyyy.ddd_yyyy.ddd
%                 'wkday'  |  1-SUN, 2-MON, ...
%                   'day'  |  yyyy.ddd
%                'tod3hr'  |  hh00 (eg. 0000, 0300, etc)
%                'tod1hr'  |  hh00 (eg. 0000, 0100, etc)
%                   '3hr'  |  yyyy.ddd.hh.mm.ss_yyyy.ddd.hh.mm.ss
%                   '1hr'  |  yyyy.ddd.hh.mm.ss_yyyy.ddd.hh.mm.ss
%
%     NEWSPAN must be one of the SPAN options (as in NOISE_STACK):
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
%     STACK2STACK(STACKDIR,NEWSPAN,'OPT1',VAL,...,'OPTN',VAL) gives access
%     to several stacking and selection parameters:
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
%     - While data i/o is .mat files by default, SAC files are okay too.
%       Note that .mat files can be read into Matlab using the LOAD
%       function and written out as SAC files using the function
%       WRITESEIZMO.  You may also convert the inputs & outputs using
%       NOISE_MAT2SAC & NOISE_SAC2MAT.  Be aware that SAC files are ignored
%       when there is a .mat file present.  Output format is based on the
%       format of the first data file read in.
%     - NEWSPAN time span must exceed the OLDSPAN time span.  For example,
%       you cannot make 3hr stacks from day stacks.  Allowed NEWSPAN values
%       given the OLDSPAN value:
%          OLDSPAN | NEWSPAN
%          --------+-------------------------------
%              all |  
%          2season | all
%          4season | 2season, all
%              mon | 4season, 2season, all
%             3mon | 4season, 2season, all
%             yrmo | mon, 3mon, 4season, 2season, all
%               wk | yrmo, mon, 3mon, 4season, 2season, all
%            wkday | all
%              day | wkday, wk, yrmo, mon, 3mon, 4season, 2season, all
%           tod3hr | all
%           tod1hr | tod3hr, all
%              3hr | tod3hr, day, wkday, wk, yrmo, mon, 3mon, 4season,
%                  |  2season, all
%              1hr | tod1hr, tod3hr, 3hr, day, wkday, wk, yrmo, mon, 3mon,
%                  |  4season, 2season, all
%
%    Header changes: SCALE (number of records in stack), DEP*
%                    Z, A, F
%
%    Examples:
%     % Creating monthly stacks from daily stacks:
%     stack2stack('outdir/zz_day_twoway','yrmo')
%
%    See also: NOISE_STACK, NOISE_SETUP, NOISE_PROCESS, NOISE_OVERVIEW
%              NOISE_STACK_ARBITRARY, NOISE_STACK_DELAZ

%     Version History:
%        Apr.  6, 2013 - initial version
%        Jan. 26, 2014 - abs path exist fix
%        May  28, 2014 - call NO_REDUNDANT_CORRELATIONS to avoid double
%                        adding due to reversed and unreversed correlations
%                        being present in a directory, remove MATIO option:
%                        now autodetects sac/mat i/o, added ZTRANSFORM
%                        option (defaults to true), added XCREVERSE option
%                        (defaults to false), bugfix: allow multiple
%                        filenames for FILENAMES option
%        May  29, 2014 - bugfix: headers/filenames now updated, bugfix:
%                        now catches error for global variable resetting
%        June  4, 2014 - also set t0 & t1 header fields
%        July  8, 2014 - edits for irlim fd i/o, set t3-4 header fields
%        July 11, 2014 - fd is converted to complex so Fisher transform
%                        works properly (FISHER was updated), iamph ok
%        July 17, 2014 - bugfix: solofun needs func handles not strings,
%                        bugfix: convert iftype to string
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 17, 2014 at 11:15 GMT

% todo:
% - overlap option
%   - default to 50% ? (current default is 100% == requires 100% overlap)
%   - weight summation by overlap?
%   - need a time function to return overlap between ranges
%   - tod?hr, wkday, mon, ?season need special care

% check nargin
error(nargchk(2,inf,nargin));
if(nargin>=3 && mod(nargin,2))
    error('seizmo:stack2stack:badInput',...
        'Unpaired option/value pair given!');
end

% directory separator
fs=filesep;

% parse/check options
opt=stack2stack_parameters(varargin{:});

% check stack directory
if(~isstring(stackdir))
    error('seizmo:stack2stack:dirNotString',...
        'STACKDIR must be a string!');
end
if(~isabspath(stackdir)); stackdir=[pwd fs stackdir]; end
if(~exist(stackdir,'dir'))
    error('seizmo:stack2stack:dirConflict',...
        ['Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end

% split stack directory
[stackpath,stackname]=fileparts(stackdir);
stackname=getwords(stackname,'_');
if(numel(stackname)~=3)
    error('seizmo:stack2stack:badInput',...
        'STACKDIR directory must be named as PAIR_OLDSPAN_XCCMP!');
end
oldspan=lower(stackname{2});

% refuse oldspan='all'
if(strcmp(oldspan,'all'))
    error('seizmo:stack2stack:badInput',...
        'STACK2STACK does not work for OLDSPAN=''ALL''!');
end

% check newspan
VALID.span_all={};
VALID.span_2season={'all'};
VALID.span_4season={'all' '2season'};
VALID.span_mon={'all' '2season' '4season'};
VALID.span_3mon={'all' '2season' '4season'};
VALID.span_yrmo={'all' '2season' '4season' '3mon' 'mon'};
VALID.span_wk={'all' '2season' '4season' '3mon' 'mon' 'yrmo'};
VALID.span_wkday={'all'};
VALID.span_day={'all' '2season' '4season' '3mon' 'mon' 'yrmo' 'wk' ...
    'wkday'};
VALID.span_tod3hr={'all'};
VALID.span_tod1hr={'all' 'tod3hr'};
VALID.span_3hr={'all' '2season' '4season' '3mon' 'mon' 'yrmo' 'wk' ...
    'wkday' 'day' 'tod3hr'};
VALID.span_1hr={'all' '2season' '4season' '3mon' 'mon' 'yrmo' 'wk' ...
    'wkday' 'day' '3hr' 'tod3hr' 'tod1hr'};
if(~ischar(newspan) || ndims(newspan)~=2 || size(newspan,1)~=1)
    error('seizmo:stack2stack:badInput',...
        'NEWSPAN must be a string!');
elseif(~ismember(lower(newspan),VALID.(['span_' oldspan])))
    error('seizmo:stack2stack:badInput',...
        'NEWSPAN=%s not allowed for OLDSPAN=%s!',newspan,oldspan);
end
newspan=lower(newspan);

% directory separator
fs=filesep;

% output directory
outdir=[stackpath fs stackname{1} '_' newspan '_' stackname{3}];

% check output directory
if(exist(outdir,'file'))
    if(~exist(outdir,'dir'))
        error('seizmo:stack2stack:dirConflict',...
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

% parallel processing setup
verbose=seizmoverbose;
%matlabpool(4); % PARALLEL

% get time-section directories
tsdirs=xdir([stackdir fs]);
tsdirs=tsdirs([tsdirs.isdir]' & ~strncmp({tsdirs.name}','.',1));
tsdirs={tsdirs.name}';

% convert input timesection names to numeric time arrays
switch oldspan
    case 'tod1hr'
        % hh00 (eg. 0000, 0100, etc)
        tmp=char(tsdirs);
        tmp=str2int(tmp(:,1:2));
        tsbgn=[ones(numel(tmp),2) tmp zeros(numel(tmp),2)];
        tsend=[ones(numel(tmp),2) tmp+1 zeros(numel(tmp),2)];
    case 'tod3hr'
        % hh00 (eg. 0000, 0300, etc)
        tmp=char(tsdirs);
        tmp=str2int(tmp(:,1:2));
        tsbgn=[ones(numel(tmp),2) tmp zeros(numel(tmp),2)];
        tsend=[ones(numel(tmp),2) tmp+3 zeros(numel(tmp),2)];
    case '1hr'
        % yyyy.ddd.hh.mm.ss_yyyy.ddd.hh.mm.ss
        tmp=char(tsdirs);
        tsbgn=[str2int(tmp(:,1:4)) str2int(tmp(:,6:8)) ...
            str2int(tmp(:,10:11)) str2int(tmp(:,13:14)) ...
            str2int(tmp(:,16:17))];
        tsend=[str2int(tmp(:,19:22)) str2int(tmp(:,24:26)) ...
            str2int(tmp(:,28:29)) str2int(tmp(:,31:32)) ...
            str2int(tmp(:,34:35))];
    case '3hr'
        % yyyy.ddd.hh.mm.ss_yyyy.ddd.hh.mm.ss
        tmp=char(tsdirs);
        tsbgn=[str2int(tmp(:,1:4)) str2int(tmp(:,6:8)) ...
            str2int(tmp(:,10:11)) str2int(tmp(:,13:14)) ...
            str2int(tmp(:,16:17))];
        tsend=[str2int(tmp(:,19:22)) str2int(tmp(:,24:26)) ...
            str2int(tmp(:,28:29)) str2int(tmp(:,31:32)) ...
            str2int(tmp(:,34:35))];
    case 'day'
        % yyyy.ddd
        tmp=char(tsdirs);
        tsbgn=[str2int(tmp(:,1:4)) str2int(tmp(:,6:8)) ...
            zeros(size(tmp,1),3)];
        tsend=tsbgn; tsend(:,2)=tsend(:,2)+1;
        tsend=fixtimes(tsend);
    case 'wkday'
        % 1-SUN, 2-MON, ...
        tmp=char(tsdirs);
        tmp=str2int(tmp(:,1));
        tsbgn=fixtimes([ones(numel(tmp),1) tmp+6 zeros(numel(tmp),3)]);
        tsend=fixtimes([ones(numel(tmp),1) tmp+7 zeros(numel(tmp),3)]);
    case 'wk'
        % yyyy.ddd_yyyy.ddd
        tmp=char(tsdirs);
        tsbgn=[str2int(tmp(:,1:4)) str2int(tmp(:,6:8)) ...
            zeros(size(tmp,1),3)];
        tsend=tsbgn; tsend(:,2)=tsend(:,2)+7;
        tsend=fixtimes(tsend);
    case 'yrmo'
        % yyyy.MM
        tmp=char(tsdirs);
        tsbgn=[str2int(tmp(:,1:4)) str2int(tmp(:,6:7)) ...
            ones(size(tmp,1),1) zeros(size(tmp,1),3)];
        tsend=tsbgn; tsend(:,2)=tsend(:,2)+1;
        tsend=fixtimes(tsend);
    case '3mon'
        % yyyy.MM_yyyy.MM
        tmp=char(tsdirs);
        tsbgn=[str2int(tmp(:,1:4)) str2int(tmp(:,6:7)) ...
            ones(size(tmp,1),1) zeros(size(tmp,1),3)];
        tsend=tsbgn; tsend(:,2)=tsend(:,2)+3;
        tsend=fixtimes(tsend);
    case 'mon'
        % MM
        tmp=str2int(char(tsdirs));
        tsbgn=[ones(numel(tmp),1) tmp ones(numel(tmp),1) ...
            zeros(numel(tmp),3)];
        tsend=tsbgn; tsend(:,2)=tsend(:,2)+1;
        tsend=fixtimes(tsend);
    case '4season'
        % spring, summer, autumn, winter
        [tsbgn,tsend]=deal(nan(numel(tsdirs),6));
        for i=1:numel(tsdirs)
            switch tsdirs{i}
                case 'spring'
                    tsbgn(i,:)=[1 1 1 0 0 0];
                    tsend(i,:)=[1 4 1 0 0 0];
                case 'summer'
                    tsbgn(i,:)=[1 4 1 0 0 0];
                    tsend(i,:)=[1 7 1 0 0 0];
                case 'autumn'
                    tsbgn(i,:)=[1 7 1 0 0 0];
                    tsend(i,:)=[1 10 1 0 0 0];
                case 'winter'
                    tsbgn(i,:)=[1 10 1 0 0 0];
                    tsend(i,:)=[2 1 1 0 0 0];
            end
        end
    case '2season'
        % summer, winter
        [tsbgn,tsend]=deal(nan(numel(tsdirs),6));
        for i=1:numel(tsdirs)
            switch tsdirs{i}
                case 'summer'
                    tsbgn(i,:)=[1 4 1 0 0 0];
                    tsend(i,:)=[1 10 1 0 0 0];
                case 'winter'
                    tsbgn(i,:)=[1 10 1 0 0 0];
                    tsend(i,:)=[2 4 1 0 0 0];
            end
        end
end

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

% determine output time spans
switch newspan
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
        spanstr=strcat(num2str(spanbgn(:,1),'%04d'),'.',...
            num2str(spanbgn(:,2),'%03d'),'.',...
            num2str(spanbgn(:,3),'%02d'),'.',...
            num2str(spanbgn(:,4),'%02d'),'.',...
            num2str(spanbgn(:,5),'%02d'),'_',...
            num2str(spanend(:,1),'%04d'),'.',...
            num2str(spanend(:,2),'%03d'),'.',...
            num2str(spanend(:,3),'%02d'),'.',...
            num2str(spanend(:,4),'%02d'),'.',...
            num2str(spanend(:,5),'%02d'));
    case '3hr'
        % 3hrs = 1/8th of a day
        spanbgn=serial2gregorian((fix(b*8):fix(e*8))/8,'doytime');
        spanend=serial2gregorian(((fix(b*8):fix(e*8))+1)/8,...
            'doytime');
        spanstr=strcat(num2str(spanbgn(:,1),'%04d'),'.',...
            num2str(spanbgn(:,2),'%03d'),'.',...
            num2str(spanbgn(:,3),'%02d'),'.',...
            num2str(spanbgn(:,4),'%02d'),'.',...
            num2str(spanbgn(:,5),'%02d'),'_',...
            num2str(spanend(:,1),'%04d'),'.',...
            num2str(spanend(:,2),'%03d'),'.',...
            num2str(spanend(:,3),'%02d'),'.',...
            num2str(spanend(:,4),'%02d'),'.',...
            num2str(spanend(:,5),'%02d'));
    case 'day'
        spanbgn=serial2gregorian(fix(b):fix(e),'doydate');
        spanend=serial2gregorian((fix(b):fix(e))+1,'doydate');
        spanstr=strcat(num2str(spanbgn(:,1),'%04d'),'.',...
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
        spanstr=strcat(num2str(spanbgn(:,1),'%04d'),'.',...
            num2str(spanbgn(:,2),'%03d'),'_',...
            num2str(spanend1(:,1),'%04d'),'.',...
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
        spanstr=strcat(num2str(spanbgn(:,1),'%04d'),'.',...
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
        spanstr=strcat(num2str(spanbgn(:,1),'%04d'),'.',...
            num2str(spanbgn(:,2),'%02d'),'_',...
            num2str(spanend1(:,1),'%04d'),'.',...
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
        spanstr=strcat(num2str(spanbgn(:,1),'%04d'),'.',...
            num2str(spanbgn(:,2),'%03d'),'.',...
            num2str(spanbgn(:,3),'%02d'),'.',...
            num2str(spanbgn(:,4),'%02d'),'.',...
            num2str(spanbgn(:,5),'%02d'),'_',...
            num2str(spanend(:,1),'%04d'),'.',...
            num2str(spanend(:,2),'%03d'),'.',...
            num2str(spanend(:,3),'%02d'),'.',...
            num2str(spanend(:,4),'%02d'),'.',...
            num2str(spanend(:,5),'%02d'));
end

% turn off checking
oldseizmocheckstate=seizmocheck_state(false);
oldcheckheaderstate=checkheader_state(false);

% for autodetecting if output is .mat or SAC files
opt.MATIO=[];

% for filetype checking
common_iftype=[];

% loop over stack time spans
%parfor s=1:size(spanbgn,1) % PARALLEL
for s=1:size(spanbgn,1) % SERIAL
    % force quietness even in parfor (which resets globals)
    seizmoverbose(false);
    
    % create variables for easy access (needed for parallel)
    sbgn=spanbgn(s,:);
    send=spanend(s,:);
    
    % get timesections ENTIRELY in stack time span
    switch newspan
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
    sdata=[];  % stack data
    snamem=[]; % stack knames
    snames=[];
    sscale=[]; % stack weights
    
    % loop over input stack dirs in this stack window
    for ts=1:nin
        % read in timesection data headers
        try
            data=load(strcat(stackdir,fs,tsdirs{in(ts)},fs,...
                'noise_records.mat'));
            data=data.noise_records;
            if(~isempty(opt.FILENAMES))
                warning('seizmo:stack2stack:unusedOption',...
                    'FILENAMES option ignored for MAT input!');
            end
            opt.MATIO_THIS_TIME=true;
            if(isempty(opt.MATIO)); opt.MATIO=true; end
        catch
            try
                data=readheader(strcat(stackdir,fs,tsdirs{in(ts)},fs,...
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
            error('seizmo:stack2stack:badInput',...
                'STACKDIR contains non-correlations!');
        elseif(numel(unique(iftype))~=1)
            error('seizmo:stack2stack:badInput',...
                'STACKDIR contains mixed correlation filetypes!');
        end
        
        % now check filetype is consistent across directories
        iftype=iftype{1};
        if(isempty(common_iftype))
            common_iftype=iftype;
        elseif(~strcmpi(common_iftype,iftype))
            error('seizmo:stack2stack:badInput',...
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
        if(verbose); print_time_left(ts,nin); end
    end
    
    % get span reference time
    switch char(newspan)
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
    sdata=changeheader(sdata,'scale',sscale,...
        'a',0,'f',timediff(sbgn,send,'utc'),...
        't0',0,'t1',timediff(sbgn,send,'utc'),...
        't3',0,'t4',timediff(sbgn,send,'utc'),...
        'z',spanref,'iztype','ia');
    
    % write out stack of stacks
    if(opt.MATIO)
        path=[outdir fs spanstr(s,:)];
        noise_records=changepath(...
            sdata,'path',path); %#ok<*NASGU>
        if(~exist(path,'dir')); mkdir(path); end
        save([path fs 'noise_records.mat'],'noise_records');
        clear noise_records;
    else
        writeseizmo(sdata,'path',path);
    end
end

% toggle checking back
seizmocheck_state(oldseizmocheckstate);
checkheader_state(oldcheckheaderstate);

% parallel processing takedown & fix verbosity
%matlabpool close; % PARALLEL
seizmoverbose(verbose);

end


function [opt]=stack2stack_parameters(varargin)
% parses/checks stack2stack pv pairs

% defaults
varargin=[{'o' 50 'ztrans' true 'q' false 'ts' [] 'te' [] 'lat' [] ...
    'lon' [] 'net' [] 'sta' [] 'str' [] 'cmp' [] 'file' [] 'xcr' false} ...
    varargin];

% require option/value pairs
if(mod(nargin,2))
    error('seizmo:stack2stack:badInput',...
        'Unpaired option/value pair given!');
elseif(~iscellstr(varargin(1:2:end)))
    error('seizmo:stack2stack:badInput',...
        'Options must be specified as strings!');
end

% get user input
for i=1:2:numel(varargin)
    switch lower(varargin{i})
        case {'overlap' 'over' 'ol' 'ov' 'o' 'olap'}
            % support for this option does not exist yet!
            if(isempty(varargin{i+1})); continue; end
            opt.OVERLAP=varargin{i+1};
        case {'xcreverse' 'xcrev' 'xcr'}
            if(isempty(varargin{i+1})); continue; end
            opt.XCREVERSE=varargin{i+1};
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
                warning('seizmo:stack2stack:badInput',...
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
            error('seizmo:stack2stack:badInput',...
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
szs=size(opt.TIMESTART);
sze=size(opt.TIMEEND);
if(~isscalar(opt.XCREVERSE) || ~islogical(opt.XCREVERSE))
    error('seizmo:stack2stack:badInput',...
        'XCREVERSE must be TRUE or FALSE!');
elseif(~isscalar(opt.OVERLAP) || ~isnumeric(opt.OVERLAP) ...
        || opt.OVERLAP>100 || opt.OVERLAP<0)
    error('seizmo:stack2stack:badInput',...
        'OVERLAP must be a scalar from 0-100%%!');
elseif(~isscalar(opt.QUIETWRITE) || ~islogical(opt.QUIETWRITE))
    error('seizmo:stack2stack:badInput',...
        'QUIETWRITE flag must be a scalar logical!');
elseif(~isempty(opt.TIMESTART) && (numel(szs)>2 || szs(1)~=1 ...
        || all(szs(2)~=[2 3 5 6]) || ~isnumeric(opt.TIMESTART) ...
        || ~isreal(opt.TIMESTART)))
    error('seizmo:stack2stack:badInput',...
        'TIMESTART must be a recognized date-time vector!');
elseif(~isempty(opt.TIMEEND) && (numel(sze)>2 || sze(1)~=1 ...
        || all(sze(2)~=[2 3 5 6]) || ~isnumeric(opt.TIMEEND) ...
        || ~isreal(opt.TIMEEND)))
    error('seizmo:stack2stack:badInput',...
        'TIMEEND must be a recognized date-time vector!');
elseif(~isempty(opt.LATRNG) && (~isnumeric(opt.LATRNG) ...
        || ~isreal(opt.LATRNG) || numel(opt.LATRNG)~=2 ...
        || size(opt.LATRNG,2)~=2 || numel(size(opt.LATRNG))~=2))
    error('seizmo:stack2stack:badInput',...
        'LATRNG must be a 2 element numeric vector as [LOW HIGH]!');
elseif(~isempty(opt.LONRNG) && (~isnumeric(opt.LONRNG) ...
        || ~isreal(opt.LONRNG) || numel(opt.LONRNG)~=2 ...
        || size(opt.LONRNG,2)~=2 || numel(size(opt.LONRNG))~=2))
    error('seizmo:stack2stack:badInput',...
        'LONRNG must be a 2 element numeric vector as [LOW HIGH]!');
elseif(~isempty(opt.NETWORKS) && (~iscellstr(opt.NETWORKS)))
    error('seizmo:stack2stack:badInput',...
        'NETWORKS must be a string list of allowed network codes!');
elseif(~isempty(opt.STATIONS) && (~iscellstr(opt.STATIONS)))
    error('seizmo:stack2stack:badInput',...
        'STATIONS must be a string list of allowed station codes!');
elseif(~isempty(opt.STREAMS) && (~iscellstr(opt.STREAMS)))
    error('seizmo:stack2stack:badInput',...
        'STREAMS must be a string list of allowed stream codes!');
elseif(~isempty(opt.COMPONENTS) && (~iscellstr(opt.COMPONENTS)))
    error('seizmo:stack2stack:badInput',...
        'COMPONENTS must be a string list of allowed component codes!');
elseif(~isempty(opt.FILENAMES) && (~iscellstr(opt.FILENAMES)))
    error('seizmo:stack2stack:badInput',...
        'FILENAMES must be a string list of allowed files!');
elseif(~isscalar(opt.ZTRANS) || ~islogical(opt.ZTRANS))
    error('seizmo:stack2stack:badInput',...
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
    error('seizmo:stack2stack:badNCFs',...
        'NCFs differ in number of points! Cannot stack!');
end
end
