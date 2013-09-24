function [lgc]=isxc(data)
%ISXC    Returns TRUE for records that are SEIZMO correlograms
%
%    Usage:    lgc=isxc(data)
%
%    Description:
%     LGC=ISXC(DATA) identifies the records in SEIZMO struct DATA that are
%     correlograms based on the header field setup as output by the
%     function CORRELATE.  LGC is equal in size to DATA and is a logical
%     array with TRUE meaning the corresponding record in DATA is a
%     correlogram.  Detection is determined by the KUSER0/1 records being
%     set to the strings 'MASTER'/'SLAVE'.
%
%    Notes:
%
%    Examples:
%     % Read in all data from a directory and exclude any correlograms:
%     d=r('*');
%     d(isxc(d))=[];
%
%    See also: CORRELATE, ROTATE_CORRELATIONS, REVERSE_CORRELATIONS,
%              HORZ_CORRELATIONS_SETS, NAME_CORRELATIONS,
%              NO_REDUNDANT_CORRELATIONS, IS_FULL_MATRIX_OF_CORRELATIONS,
%              ISSEIZMO, SEIZMOCHECK

%     Version History:
%        Feb.  5, 2013 - initial version
%        Sep. 12, 2013 - expanded see also section
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 12, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check structure
error(seizmocheck(data));

% simple checking method
% - require KUSER0=MASTER & KUSER1=SLAVE
[m,s]=getheader(data,'kuser0','kuser1');
lgc=strcmp(m,'MASTER') & strcmp(s,'SLAVE');

end
