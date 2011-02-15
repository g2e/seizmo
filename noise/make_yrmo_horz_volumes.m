function []=make_yrmo_horz_volumes(indir,startyrmo,endyrmo)
%MAKE_YRMO_HORZ_VOLUMES    Computes 1 month fk volumes for horizontals
%
%    Usage:    make_yrmo_horz_volumes(stack_dir)
%              make_yrmo_horz_volumes(stack_dir,yrmo)
%              make_yrmo_horz_volumes(stack_dir,startyrmo,endyrmo)
%
%    Description:
%     MAKE_YRMO_HORZ_VOLUMES(STACK_DIR) creates fk-based slowness response
%     volumes for an array on a month to month basis.  This only works on a
%     directory layout setup by DAYDIRS_STACKCORR.  The period range is
%     4 to 100s, the maximum slowness is 50sec/deg and the slowness
%     resolution is 1/3 sec/deg.  These are individual months, not the
%     months stacked across years of data as in MAKE_MONTHLY_HORZ_VOLUMES.
%
%     MAKE_YRMO_HORZ_VOLUMES(STACK_DIR,YRMO) explicitly sets the months to
%     compute.  Use the form YEAR.MO (eg 2006.02 without quotes).  So to
%     compute all months in 2006 use 2006.01:.01:2006.12 or use the next
%     Usage format.
%
%     MAKE_YRMO_HORZ_VOLUMES(STACK_DIR,STARTYRMO,ENDYRMO) allows setting
%     the range of months allowed using the same YRMO format above.
%
%    Notes:
%
%    Examples:
%     %
%
%    See also: MAKE_FULL_HORZ_VOLUMES, MAKE_YRMO_Z_VOLUMES, FKXCHORZVOLUME,
%              DAYDIRS_STACKCORR, MAKE_MONTHLY_HORZ_VOLUMES,
%              MAKE_DAILY_HORZ_VOLUMES

%     Version History:
%        Feb. 14, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 14, 2011 at 15:55 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% directory separator
fs=filesep;

% check stack dir
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:make_yrmo_horz_volumes:fileNotString',...
        'STACK_DIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:make_yrmo_horz_volumes:dirConflict',...
        ['STACK_DIR Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~exist([indir fs 'RR' fs 'CORR_1MONSTACK'],'dir'))
    error('seizmo:make_yrmo_horz_volumes:dirConflict',...
        ['Month Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'RR' fs 'CORR_1MONSTACK']);
end
if(~exist([indir fs 'RT' fs 'CORR_1MONSTACK'],'dir'))
    error('seizmo:make_yrmo_horz_volumes:dirConflict',...
        ['Month Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'RT' fs 'CORR_1MONSTACK']);
end
if(~exist([indir fs 'TR' fs 'CORR_1MONSTACK'],'dir'))
    error('seizmo:make_yrmo_horz_volumes:dirConflict',...
        ['Month Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'TR' fs 'CORR_1MONSTACK']);
end
if(~exist([indir fs 'TT' fs 'CORR_1MONSTACK'],'dir'))
    error('seizmo:make_yrmo_horz_volumes:dirConflict',...
        ['Month Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'TT' fs 'CORR_1MONSTACK']);
end

% get all possible yrmo
files=xdir([indir fs 'RR' fs 'CORR_1MONSTACK' fs 'STACK_*']);
names=char({files.name}');
names=str2double(cellstr(names(:,7:13)));
yrmo0=names;
files=xdir([indir fs 'RT' fs 'CORR_1MONSTACK' fs 'STACK_*']);
names=char({files.name}');
names=str2double(cellstr(names(:,7:13)));
yrmo0=[yrmo0; names];
files=xdir([indir fs 'TR' fs 'CORR_1MONSTACK' fs 'STACK_*']);
names=char({files.name}');
names=str2double(cellstr(names(:,7:13)));
yrmo0=[yrmo0; names];
files=xdir([indir fs 'TT' fs 'CORR_1MONSTACK' fs 'STACK_*']);
names=char({files.name}');
names=str2double(cellstr(names(:,7:13)));
yrmo0=unique([yrmo0; names]);

% default ranges for missing inputs
explicit=false;
if(nargin==1 || (nargin==2 && isempty(startyrmo)) ...
        || (nargin==3 && isempty(startyrmo) && isempty(endyrmo)))
    % all months available
    startyrmo=min(yrmo0);
    endyrmo=max(yrmo0);
elseif(nargin==2 || (nargin==3 && isempty(endyrmo)) ...
        || (nargin==3 && isempty(startyrmo)))
    % only months explicitly given
    yrmo=[startyrmo endyrmo];
    explicit=true;
end

% check months
if(explicit) % given months
    if(any(~isfinite(yrmo)))
        error('seizmo:make_yrmo_horz_volumes:badInput',...
            'YRMO must be a finite number of the form YEAR.MO !');
    elseif(any(round((yrmo-round(yrmo))*100)<1 ...
            | round((yrmo-round(yrmo))*100)>12))
        error('seizmo:make_yrmo_horz_volumes:badInput',...
            'YRMO must be a number of the form YEAR.MO !');
    end
else % month range
    if(~isscalar(startyrmo) || ~isscalar(endyrmo))
        error('seizmo:make_yrmo_horz_volumes:badInput',...
            'STARTYRMO & ENDYRMO must be scalar!');
    elseif(~isfinite(startyrmo) || ~isfinite(endyrmo))
        error('seizmo:make_yrmo_horz_volumes:badInput',...
            'STARTYRMO & ENDYRMO must be finite numbers!');
    elseif(round((startyrmo-round(startyrmo))*100)<1 ...
            || round((startyrmo-round(startyrmo))*100)>12)
        error('seizmo:make_yrmo_horz_volumes:badInput',...
            'STARTYRMO must be a number of the form YEAR.MO !');
    elseif(round((endyrmo-round(endyrmo))*100)<1 ...
            || round((endyrmo-round(endyrmo))*100)>12)
        error('seizmo:make_yrmo_horz_volumes:badInput',...
            'ENDYRMO must be a number of the form YEAR.MO !');
    end
end

% verbosity
verbose=seizmoverbose;
if(verbose); disp('Computing 1 month horizontal fk volumes'); end

% loop over months
for i=yrmo0(:)'
    % string version
    istr=num2str(i,'%6.2f');
    
    % skip if not in range or explicitly given
    if(explicit)
        if(~any(i==yrmo)); continue; end
    else
        if(i+eps<startyrmo || i-eps>endyrmo); continue; end
    end
    
    % detail message
    if(verbose); disp(['PROCESSING MONTH ' istr]); end
    
    % read in data
    try
        rr=readseizmo([indir fs 'RR' fs 'CORR_1MONSTACK' fs ...
            'STACK_' istr '_*_RR']);
        rt=readseizmo([indir fs 'RT' fs 'CORR_1MONSTACK' fs ...
            'STACK_' istr '_*_RT']);
        tr=readseizmo([indir fs 'TR' fs 'CORR_1MONSTACK' fs ...
            'STACK_' istr '_*_TR']);
        tt=readseizmo([indir fs 'TT' fs 'CORR_1MONSTACK' fs ...
            'STACK_' istr '_*_TT']);
    catch
        continue;
    end
    
    % find incomplete/missing stations
    [mi,si]=getheader(rr,'user0','user1');
    
    % reduce to single triangle
    rr(mi>si)=[];
    rt(mi>si)=[];
    tr(mi>si)=[];
    tt(mi>si)=[];
    
    % get fk volume
    [rvol,tvol]=fkxchorzvolume(rr,rt,tr,tt,50,301,[1/100 1/4]);
    save(['fkvol.r.' istr '.mat'],'-struct','rvol');
    save(['fkvol.t.' istr '.mat'],'-struct','tvol');
end

end
