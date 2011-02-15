function [ampl]=make_monthly_z_specampl(indir,months)
%MAKE_MONTHLY_Z_SPECAMPL    Return year-month array spectral amplitude
%
%    Usage:    ampl=make_monthly_z_specampl(stack_dir)
%              ampl=make_monthly_z_specampl(stack_dir,months)
%
%    Description:
%     AMPL=MAKE_MONTHLY_Z_SPECAMPL(STACK_DIR) returns spectral amplitudes
%     for the entire array for monthly time spans using the stacks given
%     in STACK_DIR.  Please note that these monthly time ranges span across
%     years as well to give better results for multi-year deployments.
%     STACK_DIR must have been created by DAYDIRS_STACKCORR.
%
%     AMPL=MAKE_MONTHLY_Z_SPECAMPL(STACK_DIR,MONTHS) explicitly sets the
%     months to compute.  Must be integers from 1 to 12.
%
%    Notes:
%     - AMPL is a general xy filetype.
%
%    Examples:
%     % Compare full timespan array spectral amplitudes to monthly timespan
%     % spectral amplitudes:
%     fullampl=make_full_z_specampl(my_stack_dir);
%     monthlyampl=make_monthly_z_specampl(my_stack_dir);
%     plot2([fullampl; monthlyampl]);
%
%    See also: MAKE_FULL_Z_SPECAMPL, DAYDIRS_STACKCORR,
%              MAKE_DAILY_Z_SPECAMPL, MAKE_YRMO_Z_SPECAMPL

%     Version History:
%        Feb. 13, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 13, 2011 at 15:55 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% directory separator
fs=filesep;

% check stack dir
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:make_monthly_z_specampl:fileNotString',...
        'STACK_DIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:make_monthly_z_specampl:dirConflict',...
        ['STACK_DIR Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~exist([indir fs 'ZZ' fs 'CORR_MONSTACK'],'dir'))
    error('seizmo:make_monthly_z_specampl:dirConflict',...
        ['Month Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'ZZ' fs 'CORR_MONSTACK']);
end

% default/check months
if(nargin==1 || isempty(months)); months=1:12; end
if(~isreal(months) || any(months(:)~=fix(months(:))) ...
        || any(months(:)<1 | months(:)>12))
    error('seizmo:make_monthly_z_specampl:badMonths',...
        'MONTHS must be an array of integers within 1 & 12!');
end

% loop over months
cnt=1;
month={'JAN' 'FEB' 'MAR' 'APR' 'MAY' 'JUN' ...
    'JUL' 'AUG' 'SEP' 'OCT' 'NOV' 'DEC'};
for i=months
    % read in data
    try
        data=readseizmo([indir fs 'ZZ' fs 'CORR_MONSTACK' fs ...
            'STACK_' num2str(i,'%02d') '_*_ZZ']);
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
        'kname',nan,'scale',numel(data),'kstnm',month{i});
    
    % change name
    ampl(cnt,1)=changename(ampl(cnt,1),'name',...
        ['ARRAY_SPECAMPL_' num2str(i,'%02d') '_ZZ']);
    ampl(cnt,1)=changepath(ampl(cnt,1),'path','.');
    
    % increment
    cnt=cnt+1;
end

end
