function [ampl]=make_yrmo_z_specampl(indir,startyrmo,endyrmo)
%MAKE_YRMO_Z_SPECAMPL    Return monthly array spectral amplitude
%
%    Usage:    ampl=make_yrmo_z_specampl(stack_dir)
%              ampl=make_yrmo_z_specampl(stack_dir,yrmo)
%              ampl=make_yrmo_z_specampl(stack_dir,startyrmo,endyrmo)
%
%    Description:
%     AMPL=MAKE_YRMO_Z_SPECAMPL(STACK_DIR) returns spectral amplitudes
%     for the entire array for monthly time spans using the stacks given
%     in STACK_DIR.  These are individual months, not the months stacked
%     across years of data.  STACK_DIR must have been created by
%     DAYDIRS_STACKCORR.
%
%     AMPL=MAKE_YRMO_Z_SPECAMPL(STACK_DIR,YRMO) explicitly sets the
%     months to compute.  Use the form YEAR.MO (eg 2006.02 without quotes).
%     So to compute all months in 2006 use: 2006.01:.01:2006.12 or use the
%     next Usage format.
%
%     AMPL=MAKE_YRMO_Z_SPECAMPL(STACK_DIR,STARTYRMO,ENDYRMO) allows setting
%     the range of months allowed using the same YRMO format above.
%
%    Notes:
%     - AMPL is a general xy filetype.
%
%    Examples:
%     % 
%
%    See also: MAKE_FULL_Z_SPECAMPL, DAYDIRS_STACKCORR,
%              MAKE_MONTHLY_Z_SPECAMPL, MAKE_DAILY_Z_SPECAMPL

%     Version History:
%        Feb. 13, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 13, 2011 at 15:55 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% directory separator
fs=filesep;

% check stack dir
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:make_yrmo_z_specampl:fileNotString',...
        'STACK_DIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:make_yrmo_z_specampl:dirConflict',...
        ['STACK_DIR Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~exist([indir fs 'ZZ' fs 'CORR_1MONSTACK'],'dir'))
    error('seizmo:make_yrmo_z_specampl:dirConflict',...
        ['Month Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'ZZ' fs 'CORR_1MONSTACK']);
end

% get all possible yrmo
files=xdir([indir fs 'ZZ' fs 'CORR_1MONSTACK' fs 'STACK_*']);
names=char({files.name}');
names=str2double(cellstr(names(:,7:13)));
yrmo0=unique(names);

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
        error('seizmo:make_yrmo_z_specampl:badInput',...
            'YRMO must be a finite number of the form YEAR.MO !');
    elseif(any(round((yrmo-round(yrmo))*100)<1 ...
            | round((yrmo-round(yrmo))*100)>12))
        error('seizmo:make_yrmo_z_specampl:badInput',...
            'YRMO must be a number of the form YEAR.MO !');
    end
else % month range
    if(~isscalar(startyrmo) || ~isscalar(endyrmo))
        error('seizmo:make_yrmo_z_specampl:badInput',...
            'STARTYRMO & ENDYRMO must be scalar!');
    elseif(~isfinite(startyrmo) || ~isfinite(endyrmo))
        error('seizmo:make_yrmo_z_specampl:badInput',...
            'STARTYRMO & ENDYRMO must be finite numbers!');
    elseif(round((startyrmo-round(startyrmo))*100)<1 ...
            || round((startyrmo-round(startyrmo))*100)>12)
        error('seizmo:make_yrmo_z_specampl:badInput',...
            'STARTYRMO must be a number of the form YEAR.MO !');
    elseif(round((endyrmo-round(endyrmo))*100)<1 ...
            || round((endyrmo-round(endyrmo))*100)>12)
        error('seizmo:make_yrmo_z_specampl:badInput',...
            'ENDYRMO must be a number of the form YEAR.MO !');
    end
end

% verbosity (turn it off for the loop)
verbose=seizmoverbose(false);
if(verbose); disp('Computing monthly array amplitude spectra'); end

% loop over months
cnt=1;
for i=yrmo0(:)'
    % string version
    si=num2str(i);
    
    % skip if not in range or explicitly given
    if(explicit)
        if(~any(i==yrmo)); continue; end
    else
        if(i+eps<startyrmo || i-eps>endyrmo); continue; end
    end
    
    % detail message
    if(verbose); disp(['PROCESSING MONTH ' si]); end
    
    % read in data
    try
        data=readseizmo([indir fs 'ZZ' fs 'CORR_1MONSTACK' fs ...
            'STACK_' si '_*_ZZ']);
    catch
        continue;
    end
    
    % get number of days in each stack
    days=getheader(data,'scale');
    days(isnan(days))=1;
    
    % weighted average
    ampl(cnt,1)=divide(stack(keepam(dft(multiply(data,days)))),sum(days));
    
    % clean headers somewhat
    ampl(cnt,1)=changeheader(ampl(cnt,1),'st',nan,'ev',nan,'cmp',nan,...
        'kname',nan,'scale',numel(data),'kstnm',si);
    
    % change name
    ampl(cnt,1)=changename(ampl(cnt,1),'name',...
        ['ARRAY_SPECAMPL_' si '_ZZ']);
    ampl(cnt,1)=changepath(ampl(cnt,1),'path','.');
    
    % increment
    cnt=cnt+1;
end

% return verbosity to normal
seizmoverbose(verbose);

end
