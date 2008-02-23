function [data]=fixdelta(data)
%FIXDELTA    Fix sample spacing for SAClab data records
%
%    Description:  Fixes the sample spacing of data records to be the
%     decimal equivalent of a fraction of 2 small integers.  This is
%     particularly useful for upgrading the sample spacing from single to
%     double precision.
%
%    Usage: data=fixdelta(data)
%
%    Examples:
%     Double the precision of records and fix the delta intervals
%      data=fixdelta(doubleit(data))
%
%    See also: doubleit, chkhdr

[n,d]=rat(gh(data,'delta'));
data=ch(data,'delta',n./d);

end

