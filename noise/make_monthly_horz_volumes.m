function []=make_monthly_horz_volumes(indir,months)
%MAKE_MONTHLY_HORZ_VOLUMES    Computes fk volumes for horiz. monthly stacks
%
%    Usage:    make_monthly_horz_volumes(stack_dir)
%              make_monthly_horz_volumes(stack_dir,months)
%
%    Description: MAKE_MONTHLY_HORZ_VOLUMES(STACK_DIR) creates fk-based
%     slowness response volumes for an array on a month to month basis.
%     This only works on a directory layout setup by DAYDIRS_STACKCORR.
%     The period range is 4 to 100s, the maximum slowness is 50sec/deg and
%     the slowness resolution is 1/3 sec/deg.  This will process all 12
%     months of the year.
%
%     MAKE_MONTHLY_HORZ_VOLUMES(STACK_DIR,MONTHS) explicitly sets the
%     months to process.  This is useful for doing this parallel-ish.
%
%    Notes:
%
%    Examples:
%
%    See also: MAKE_FULL_HORZ_VOLUMES, MAKE_MONTHLY_Z_VOLUMES,
%              FKXCHORZVOLUME, DAYDIRS_STACKCORR

%     Version History:
%        June 24, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 24, 2010 at 15:55 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% directory separator
fs=filesep;

% check stack dir
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:make_monthly_horz_volumes:fileNotString',...
        'STACK_DIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:make_monthly_horz_volumes:dirConflict',...
        ['STACK_DIR Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~exist([indir fs 'RR' fs 'CORR_MONSTACK'],'dir'))
    error('seizmo:make_monthly_horz_volumes:dirConflict',...
        ['Month Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'RR' fs 'CORR_MONSTACK']);
end
if(~exist([indir fs 'RT' fs 'CORR_MONSTACK'],'dir'))
    error('seizmo:make_monthly_horz_volumes:dirConflict',...
        ['Month Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'RT' fs 'CORR_MONSTACK']);
end
if(~exist([indir fs 'TR' fs 'CORR_MONSTACK'],'dir'))
    error('seizmo:make_monthly_horz_volumes:dirConflict',...
        ['Month Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'TR' fs 'CORR_MONSTACK']);
end
if(~exist([indir fs 'TT' fs 'CORR_MONSTACK'],'dir'))
    error('seizmo:make_monthly_horz_volumes:dirConflict',...
        ['Month Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'TT' fs 'CORR_MONSTACK']);
end

% default/check months
if(nargin==1 || isempty(months)); months=1:12; end
if(~isreal(months) || any(months(:)~=fix(months(:))) ...
        || any(months(:)<1 || months(:)>12))
    error('seizmo:make_monthly_horz_volumes:badMonths',...
        'MONTHS must be an array of integers within 1 & 12!');
end

% loop over months
for i=months
    % read in data
    rr=readseizmo([indir fs 'RR' fs 'CORR_MONSTACK' fs ...
        'STACK_' num2str(i,'%02d') '_*']);
    rt=readseizmo([indir fs 'RT' fs 'CORR_MONSTACK' fs ...
        'STACK_' num2str(i,'%02d') '_*']);
    tr=readseizmo([indir fs 'TR' fs 'CORR_MONSTACK' fs ...
        'STACK_' num2str(i,'%02d') '_*']);
    tt=readseizmo([indir fs 'TT' fs 'CORR_MONSTACK' fs ...
        'STACK_' num2str(i,'%02d') '_*']);
    
    % find incomplete/missing stations
    [mi,si]=getheader(rr,'user0','user1');
    
    % reduce to single triangle
    rr(mi>si)=[];
    rt(mi>si)=[];
    tr(mi>si)=[];
    tt(mi>si)=[];
    
    % get fk volume
    [rvol,tvol]=fkxchorzvolume(rr,rt,tr,tt,50,301,[1/100 1/4]);
    save(['svol.r.' num2str(i,'%02d') '.mat'],'-struct','rvol');
    save(['svol.t.' num2str(i,'%02d') '.mat'],'-struct','tvol');
end

end
