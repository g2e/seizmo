function []=make_full_z_volumes(indir)
%MAKE_FULL_Z_VOLUMES    Computes fk volumes for vertical full stacks
%
%    Usage:    make_full_z_volumes(stack_dir)
%
%    Description: MAKE_FULL_Z_VOLUMES(STACK_DIR) creates a fk-based
%     slowness response volume for an array.  This only works on a
%     directory layout setup by DAYDIRS_STACKCORR.  The period range is 4
%     to 100s, the maximum slowness is 50sec/deg and the slowness
%     resolution is 1/3 sec/deg.  This will process only the full stack.
%
%    Notes:
%
%    Examples:
%
%    See also: MAKE_MONTHLY_Z_VOLUMES, MAKE_FULL_HORZ_VOLUMES,
%              FKXCVOLUME, DAYDIRS_STACKCORR

%     Version History:
%        June 24, 2010 - initial version
%        Oct. 10, 2010 - svol to fkvol
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 10, 2010 at 15:55 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% directory separator
fs=filesep;

% check stack dir
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:make_full_z_volumes:fileNotString',...
        'STACK_DIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:make_full_z_volumes:dirConflict',...
        ['STACK_DIR Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~exist([indir fs 'ZZ' fs 'CORR_FULLSTACK'],'dir'))
    error('seizmo:make_full_z_volumes:dirConflict',...
        ['Full Stack Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],...
        [indir fs 'ZZ' fs 'CORR_FULLSTACK']);
end

% read in data
data=readseizmo([indir fs 'ZZ' fs 'CORR_FULLSTACK' fs 'STACK_FULL_*']);

% find incomplete/missing stations
[mi,si]=getheader(data,'user0','user1');

% reduce to single triangle
data(mi>si)=[];

% get fk volume
vol=fkxcvolume(data,50,301,[1/100 1/4]);
save('fkvol.z.fullstack.mat','-struct','vol');

end
