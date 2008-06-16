function [data]=sub(data,constant,cmp)
%SUB    Subtract a constant from SAClab data records
%
%    Description: SUB(DATA,CONSTANT) subtracts a constant from SAClab data
%     records.  For multi-component files, this operation is performed on
%     every component (this includes spectral files).
%
%     SUB(DATA,CONSTANT,CMP) allows for operations on just components in
%     the list CMP.  See the examples section for a usage case.
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Usage: data=sub(data,constant)
%           data=sub(data,constant,cmp_list)
%
%    Notes:
%     - a scalar constant applies the value to all records
%     - a vector of constants (length must equal the number of records)
%       allows applying different values to each record
%     - cmp_list gives the dependent component(s) to work on (default=all)
%
%    Examples:
%     Do a Hilbert transform by converting records to the frequency 
%     domain, subtracting pi/2 from the phase (component 2 in spectral
%     records), and converting back to the time domain:
%      data=idft(sub(dft(data),pi/2,2))
%
%    See also: add, mul, divide

%     Version History:
%        ????????????? - Initial Version
%        June 11, 2008 - Cleaned up documentation and added example
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 11, 2008 at 24:00 GMT

% check nargin
error(nargchk(2,3,nargin))

% no constant case
if(isempty(constant)); return; end

% default component
if(nargin==2 || isempty(cmp)); cmp=':'; end

% send to add
data=add(data,-constant,cmp);

end
