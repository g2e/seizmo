function []=make_monthly_z_volumes(indir,months)
%MAKE_MONTHLY_Z_VOLUMES    Computes year-month fk volumes for verticals
%
%    Usage:    make_monthly_z_volumes(stack_dir)
%              make_monthly_z_volumes(stack_dir,months)
%
%    Description:
%     MAKE_MONTHLY_Z_VOLUMES(STACK_DIR) creates fk-based slowness response
%     volumes for an array on a month to month basis independent of year
%     (eg stacked across years).  This only works on a directory layout
%     setup by DAYDIRS_STACKCORR.  The period range is 4 to 100s, the
%     maximum slowness is 50sec/deg and the slowness resolution is 1/3
%     sec/deg.
%
%     MAKE_MONTHLY_Z_VOLUMES(STACK_DIR,MONTHS) explicitly sets the months
%     to process.  This is useful for doing this parallel-ish.
%
%    Notes:
%
%    Examples:
%     % 
%
%    See also: MAKE_FULL_Z_VOLUMES, MAKE_MONTHLY_HORZ_VOLUMES, FKXCVOLUME,
%              DAYDIRS_STACKCORR, MAKE_YRMO_Z_VOLUMES, MAKE_DAILY_Z_VOLUMES

%     Version History:
%        June 24, 2010 - initial version
%        Oct. 10, 2010 - svol to fkvol
%        Feb. 14, 2011 - doc update, skip if no data for a month, verbosity
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 15:55 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% directory separator
fs=filesep;

% check stack dir
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:make_monthly_z_volumes:fileNotString',...
        'STACK_DIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:make_monthly_z_volumes:dirConflict',...
        ['STACK_DIR Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~exist([indir fs 'ZZ' fs 'CORR_MONSTACK'],'dir'))
    error('seizmo:make_monthly_z_volumes:dirConflict',...
        ['Month Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'ZZ' fs 'CORR_MONSTACK']);
end

% default/check months
if(nargin==1 || isempty(months)); months=1:12; end
if(~isreal(months) || any(months(:)~=fix(months(:))) ...
        || any(months(:)<1 | months(:)>12))
    error('seizmo:make_monthly_z_volumes:badMonths',...
        'MONTHS must be an array of integers within 1 & 12!');
end

% verbosity
verbose=seizmoverbose;
if(verbose); disp('Computing monthly fk volumes'); end

% loop over months
for i=months
    % detail message
    if(verbose); disp(['PROCESSING MONTH ' num2str(i,'%02d')]); end
    
    % read in data
    try
        data=readseizmo([indir fs 'ZZ' fs 'CORR_MONSTACK' fs ...
            'STACK_' num2str(i,'%02d') '_*_ZZ']);
    catch
        continue;
    end
    
    % find incomplete/missing stations
    [mi,si]=getheader(data,'user0','user1');
    
    % reduce to single triangle
    data(mi>si)=[];
    
    % get fk volume
    vol=fkxcvolume(data,50,301,[1/100 1/4]);
    save(['fkvol.z.' num2str(i,'%02d') '.mat'],'-struct','vol');
end

end
