function [bytes]=seizmosize(data)
%SEIZMOSIZE    Returns header-estimated disksize of SEIZMO records in bytes
%
%    Usage:    bytes=seizmosize(data)
%
%    Description: SEIZMOSIZE(DATA) returns the expected on-disk size in
%     bytes of each record in DATA using only the header info.  This is
%     mainly to detect files on disk that are inconsistent in size from
%     that expected given their header info.
%
%    Notes:
%
%    Examples:
%     Calculate expected filesizes for SAC files in current directory:
%      seizmosize(readheader('*.SAC'))
%
%    See also: readheader, readdata

%     Version History:
%        Feb. 29, 2008 - initial version (was a subfunction called sacsize)
%        Mar.  4, 2008 - doc update, use LGCCHK
%        June 12, 2008 - doc update
%        Sep. 14, 2008 - doc update
%        Sep. 25, 2008 - GET_N_CHECK adds dataless support
%        Sep. 26, 2008 - VINFO cleans up/reduces code
%        Oct.  7, 2008 - drop GET_N_CHECK (keep dataless support)
%        Oct. 26, 2008 - CHKHDR added (true dataless support), doc update
%        Nov. 13, 2008 - renamed from SEISSIZE to SEIZSIZE
%        Nov. 15, 2008 - update for new naming scheme (now SEIZMOSIZE)
%        Apr. 23, 2009 - fix nargchk for octave
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 21:15 GMT

% todo:

% check number of inputs
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% headers setup (check struct too)
[h,vi]=versioninfo(data);
hdata=[h(vi).data];

% turn off struct checking
oldstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check header
data=checkheader(data);

% header info
ncmp=getncmp(data);
npts=getheader(data,'npts');
iftype=getenumdesc(data,'iftype');
leven=getlgc(data,'leven');

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

% toggle struct checking back
set_seizmocheck_state(oldstate);

end
