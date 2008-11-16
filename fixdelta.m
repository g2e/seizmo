function [data]=fixdelta(data)
%FIXDELTA    Fix sample spacing for SEIZMO data records
%
%    Description:  FIXDELTA(DATA) fixes the sample spacing (DELTA header
%     field) of records in DATA to be the decimal equivalent of a fraction
%     of 2 small integers.  This is particularly useful for upgrading the
%     sample spacing from single to double precision as it extends the
%     precision.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Header changes: DELTA
%
%    Usage:    data=fixdelta(data)
%
%    Examples:
%     Force the double precision and update the delta field:
%      data=fixdelta(recast(data,'double'));
%
%    See also: rat, recast

%     Version History:
%        Feb. 21, 2008 - initial version
%        Feb. 23, 2008 - minor doc update
%        Feb. 28, 2008 - minor doc update
%        Mar.  4, 2008 - minor doc update
%        Oct.  8, 2008 - doc update, add history
%        Nov. 16, 2008 - update for new name schema, doc and history update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 16, 2008 at 06:10 GMT

% todo:

[n,d]=rat(getheader(data,'delta'));
data=changeheader(data,'delta',n./d);

end
