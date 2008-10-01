function [bytes]=seissize(data)
%SEISSIZE    Returns header estimated size of SAClab datafiles in bytes
%
%    Description: SEISSIZE(DATA) returns the expected size (in bytes) using
%     the header info of each datafile that would be written from DATA. 
%     This is mainly used to rapidly detect files on disk that are 
%     inconsistent in size from what their header info indicates.
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Header changes: N/A
%
%    Usage:  bytes=seissize(data)
%
%    Examples:
%     Calculate expected filesizes for SAC files in current directory:
%      seissize(rh('*.SAC'))
%
%    See also: rh, rdata

%     Version History:
%        Feb. 29, 2008 - initial version (from subfunction sacsize)
%        Mar.  4, 2008 - doc update, use LGCCHK
%        June 12, 2008 - doc update
%        Sep. 14, 2008 - doc update
%        Sep. 25, 2008 - GET_N_CHECK adds dataless support
%        Sep. 26, 2008 - VINFO cleans up/reduces code
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 26, 2008 at 19:45 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data))

% pull necessary header info
[ncmp,npts,iftype,leven]=get_n_check(data);

% headers setup
[h,vi]=vinfo(data);
hdata=[h(vi).data];

% filetype
count=strcmpi(iftype,'Time Series File')...
    +strcmpi(iftype,'General X vs Y file')...
    +strcmpi(iftype,'General XYZ (3-D) file')...
    +2*strcmpi(iftype,'Spectral File-Real/Imag')...
    +2*strcmpi(iftype,'Spectral File-Ampl/Phase');

% multi-component
count=ncmp.*count;

% uneven sampling
count=count+strcmpi(leven,'false');

% final tally
bytes=[hdata.startbyte].'+count.*npts.*[hdata.bytesize].';

end
