function [data]=fixdelta(data)
%FIXDELTA    Fix sample spacing for SAClab data records
%
%    Description:  FIXDELTA(DATA) fixes the sample spacing (DELTA header
%     field) of records in DATA to be the decimal equivalent of a fraction
%     of 2 small integers.  This is particularly useful for upgrading the
%     sample spacing from single to double precision as it extends the
%     precision.
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Header changes: DELTA
%
%    Usage: data=fixdelta(data)
%
%    Examples:
%     Double the precision of records and fix the delta intervals
%      data=fixdelta(classit(data,'double'));
%
%    See also: doubleit, chkhdr

%     Version History:
%        Feb. 21, 2008 - initial version
%        Feb. 23, 2008 - uses GLGC now
%        Feb. 28, 2008 - uses SEISCHK now
%        Mar.  4, 2008 - doc update, fixed LEVEN bug, uses LGCCHK now
%        Oct.  8, 2008 - doc update, add history
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct.  8, 2008 at 06:50 GMT

% todo:

[n,d]=rat(gh(data,'delta'));
data=ch(data,'delta',n./d);

end
