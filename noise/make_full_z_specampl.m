function [ampl]=make_full_z_specampl(indir)
%MAKE_FULL_Z_SPECAMPL    Return full-time array spectral amplitude
%
%    Usage:    ampl=make_full_z_specampl(stack_dir)
%
%    Description:
%     AMPL=MAKE_FULL_Z_SPECAMPL(STACK_DIR) returns the spectral amplitude
%     for the entire array for the entire time span using the stacks given
%     in STACK_DIR.  STACK_DIR must have been created by DAYDIRS_STACKCORR.
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
%    See also: MAKE_MONTHLY_Z_SPECAMPL, DAYDIRS_STACKCORR,
%              MAKE_DAILY_Z_SPECAMPL, MAKE_YRMO_Z_SPECAMPL

%     Version History:
%        Feb. 13, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 13, 2011 at 15:55 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% directory separator
fs=filesep;

% check stack dir
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:make_full_z_specampl:fileNotString',...
        'STACK_DIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:make_full_z_specampl:dirConflict',...
        ['STACK_DIR Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~exist([indir fs 'ZZ' fs 'CORR_FULLSTACK'],'dir'))
    error('seizmo:make_full_z_specampl:dirConflict',...
        ['Full Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'ZZ' fs 'CORR_FULLSTACK']);
end

% read in data
data=readseizmo([indir fs 'ZZ' fs 'CORR_FULLSTACK' fs 'STACK_FULL_*_ZZ']);

% get number of days in each stack
days=getheader(data,'scale');
days(isnan(days))=1;

% weighted average
ampl=divide(stack(keepam(dft(multiply(data,days)))),sum(days));

% clean headers somewhat
ampl=changeheader(ampl,'st',nan,'ev',nan,'cmp',nan,'kname',nan,...
    'scale',numel(data),'kstnm','FULL');

% change name
ampl=changename(ampl,'name','ARRAY_SPECAMPL_FULL_ZZ');
ampl=changepath(ampl,'path','.');

end
