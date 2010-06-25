function []=make_full_horz_volumes(indir)
%MAKE_FULL_HORZ_VOLUMES    Computes fk volumes for horizontal full stacks
%
%    Usage:    make_full_horz_volumes(stack_dir)
%
%    Description: MAKE_FULL_HORZ_VOLUMES(STACK_DIR) creates a fk-based
%     slowness response volume for an array.  This only works on a
%     directory layout setup by DAYDIRS_STACKCORR.  The period range is 4
%     to 100s, the maximum slowness is 50sec/deg and the slowness
%     resolution is 1/3 sec/deg.  This will process only the full stack.
%
%    Notes:
%
%    Examples:
%
%    See also: MAKE_MONTHLY_HORZ_VOLUMES, MAKE_FULL_Z_VOLUMES,
%              FKXCHORZVOLUME, DAYDIRS_STACKCORR

%     Version History:
%        June 24, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 24, 2010 at 15:55 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% directory separator
fs=filesep;

% check stack dir
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:make_full_horz_volumes:fileNotString',...
        'STACK_DIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:make_full_horz_volumes:dirConflict',...
        ['STACK_DIR Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~exist([indir fs 'RR' fs 'CORR_FULLSTACK'],'dir'))
    error('seizmo:make_full_horz_volumes:dirConflict',...
        ['Full Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'RR' fs 'CORR_FULLSTACK']);
end
if(~exist([indir fs 'RT' fs 'CORR_FULLSTACK'],'dir'))
    error('seizmo:make_full_horz_volumes:dirConflict',...
        ['Full Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'RT' fs 'CORR_FULLSTACK']);
end
if(~exist([indir fs 'TR' fs 'CORR_FULLSTACK'],'dir'))
    error('seizmo:make_full_horz_volumes:dirConflict',...
        ['Full Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'TR' fs 'CORR_FULLSTACK']);
end
if(~exist([indir fs 'TT' fs 'CORR_FULLSTACK'],'dir'))
    error('seizmo:make_full_horz_volumes:dirConflict',...
        ['Full Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'TT' fs 'CORR_FULLSTACK']);
end

% read in data
rr=readseizmo([indir fs 'RR' fs 'CORR_FULLSTACK' fs 'STACK_FULL_*']);
rt=readseizmo([indir fs 'RT' fs 'CORR_FULLSTACK' fs 'STACK_FULL_*']);
tr=readseizmo([indir fs 'TR' fs 'CORR_FULLSTACK' fs 'STACK_FULL_*']);
tt=readseizmo([indir fs 'TT' fs 'CORR_FULLSTACK' fs 'STACK_FULL_*']);

% find incomplete/missing stations
[mi,si]=getheader(rr,'user0','user1');

% reduce to single triangle
rr(mi>si)=[];
rt(mi>si)=[];
tr(mi>si)=[];
tt(mi>si)=[];

% get fk volume
[rvol,tvol]=fkxchorzvolume(rr,rt,tr,tt,50,301,[1/100 1/4]);
save('svol.r.fullstack.mat','-struct','rvol');
save('svol.t.fullstack.mat','-struct','tvol');

end
