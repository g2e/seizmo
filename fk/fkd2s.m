function [fk]=fkd2s(fk)
%FKD2S    Convert fk response from double to single precision
%
%    Usage:    fkstruct=fkd2s(fkstruct)
%
%    Description: FKSTRUCT=FKD2S(FKSTRUCT) converts the response in
%     FK struct FKSTRUCT (must be compatible with FKMAP/FKVOLUME/FK4D
%     output) from double precision to single precision.  This reduces
%     memory consumption significantly and is now the default for output
%     from FK functions.
%
%    Notes:
%
%    Examples:
%     Read in a fk volume made in double precision, convert to single
%     and write back out:
%      vol=load('myfkvol.mat');
%      vol=fkd2s(vol);
%      save('myfkvol.mat','-struct','vol');
%
%    See also: FKMAP, FKVOLUME, FK4D, FKXCVOLUME, FKHORZXCVOLUME

%     Version History:
%        June 11, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 11, 2010 at 14:55 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check fk struct
error(chkfkstruct(fk));

% convert double to single
for i=1:numel(fk)
    fk(i).response=single(fk(i).response);
end

end
