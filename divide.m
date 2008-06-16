function [data]=divide(data,constant,cmp)
%DIVIDE    Divide SAClab data records by a constant
%
%    Description: DIVIDE(DATA,CONSTANT) divides SAClab data records by a
%     constant.  For multi-component files, this operation is performed on
%     every component (this includes spectral files).
%
%     DIVIDE(DATA,CONSTANT,CMP) allows for operations on just components in
%     the list CMP.  See the examples section for a usage case.
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Usage: data=divide(data,constant)
%           data=divide(data,constant,cmp_list)
%
%    Notes:
%     - a scalar constant applies the value to all records
%     - a vector of constants (length must equal the number of records)
%       allows applying different values to each record
%     - cmp_list gives the dependent component(s) to work on (default=all)
%
%    Examples:
%     Alter the amplitudes of amplitude-phase spectral records without
%     affecting the phase component by dividing only the first component:
%      data=divide(data,32,1)
%
%    See also: mul, add, sub

%     Version History:
%        ????????????? - Initial Version
%        June 12, 2008 - Cleaned up documentation and added example
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2008 at 00:15 GMT

% check nargin
error(nargchk(2,3,nargin))

% no constant case
if(isempty(constant)); return; end

% default component
if(nargin==2 || isempty(cmp)); cmp=':'; end

% send to mult
data=mult(data,1./constant,cmp);

end
