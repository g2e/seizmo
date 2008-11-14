function [bytes]=seizsize(data)
%SEIZSIZE    Returns header-estimated disksize of SAClab datafiles in bytes
%
%    Description: SEIZSIZE(DATA) returns the expected on-disk size in bytes
%     of each record in DATA using only the header info.  This is mainly
%     to rapidly detect files on disk that are inconsistent in size from
%     that expected given their header info.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Usage:    bytes=seizsize(data)
%
%    Examples:
%     Calculate expected filesizes for SAC files in current directory:
%      seizsize(rh('*.SAC'))
%
%    See also: readhdr, readdata

%     Version History:
%        Feb. 29, 2008 - initial version (from subfunction sacsize)
%        Mar.  4, 2008 - doc update, use LGCCHK
%        June 12, 2008 - doc update
%        Sep. 14, 2008 - doc update
%        Sep. 25, 2008 - GET_N_CHECK adds dataless support
%        Sep. 26, 2008 - VINFO cleans up/reduces code
%        Oct.  7, 2008 - drop GET_N_CHECK (keep dataless support)
%        Oct. 26, 2008 - CHKHDR added (true dataless support), doc update
%        Nov. 13, 2008 - renamed from SEISSIZE to SEIZSIZE
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 13, 2008 at 04:30 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin))

% headers setup (check struct too)
[h,vi]=versioninfo(data);
hdata=[h(vi).data];

% turn off struct checking
oldstate=get_seizchk_state;
set_seizchk_state(false);

% check header
data=chkhdr(data);

% header info
ncmp=getncmp(data);
npts=gethdr(data,'npts');
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
set_seizchk_state(oldstate);

end
