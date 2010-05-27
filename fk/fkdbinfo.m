function [mindb,meddb,maxdb]=fkdbinfo(fk)
%FKDBINFO    Returns the min/median/max dB for a FK struct
%
%    Usage:    [mindb,meddb,maxdb]=fkdbinfo(fk)
%
%    Description: [MINDB,MEDDB,MAXDB]=FKDBINFO(FK) returns the limits and
%     median of the response(s) in FK in decibels.  This is useful for
%     quick identification of elements with strong plane wave coherency.
%
%    Notes:
%
%    Examples:
%     Analyze a 4D fk dataset:
%      s4d=fk4d(data,[],[],50,201,[1/100 1/4]);
%      [mindb,meddb,maxdb]=fkdbinfo(s4d);
%      figure; plot(mindb);
%      hold on;
%      plot(maxdb);
%      plot(meddb);
%
%    See also: FKMAP, FKVOLUME, FK4D

%     Version History:
%        May  26, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  26, 2010 at 09:15 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check fk struct
error(chkfkstruct(fk));

% loop over every element
mindb=nan(size(fk));
meddb=mindb;
maxdb=mindb;
for i=1:numel(fk)
    mindb(i)=min(fk(i).response(:))+fk(i).normdb;
    meddb(i)=median(fk(i).response(:))+fk(i).normdb;
    maxdb(i)=max(fk(i).response(:))+fk(i).normdb;
end

end
