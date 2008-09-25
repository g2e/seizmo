function [ncmp,npts,iftype,leven]=get_n_check(data)
%GET_N_CHECK    Gets and checks some basic SAClab header fields
%
%    Description: GET_N_CHECK(DATA) returns the values associated with
%     several SAClab header fields and also does some common checks along
%     the way.  This is an internal function to reduce rampant code 
%     repetition.
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Input/Output requirements: DATA must be a valid SAClab structure
%
%    Header changes: N/A
%
%    Usage:    [ncmp,npts,iftype,leven]=get_n_check(data)
%
%    Examples:
%     NONE
%
%    See also:

%     Version History:
%        Sep. 25, 2008 - initial version, adds dataless support
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 25, 2008 at 06:05 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data))

% header info
iftype=genumdesc(data,'iftype');
warning('off','SAClab:gh:fieldInvalid')
[npts,ncmp]=gh(data,'npts','ncmp');
warning('on','SAClab:gh:fieldInvalid')
leven=glgc(data,'leven');
error(lgcchk('leven',leven(npts>1)))

% clean up and check ncmp
ncmp(isnan(ncmp))=1;
if(any(ncmp<1 | fix(ncmp)~=ncmp))
    error('SAClab:get_n_check:badNumCmp',...
        'Field NCMP must be a positive integer!')
end

end