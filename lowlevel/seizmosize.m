function [bytes,hbytes,dbytes]=seizmosize(data)
%SEIZMOSIZE    Returns header-estimated disksize of SEIZMO records in bytes
%
%    Usage:    bytes=seizmosize(data)
%              [bytes,hbytes,dbytes]=seizmosize(data)
%
%    Description:
%     BYTES=SEIZMOSIZE(DATA) returns the expected on-disk size in bytes for
%     each record in DATA using only the header info.  This is mainly to
%     detect files on disk that are inconsistent in size from that expected
%     given their header info.
%
%    [BYTES,HBYTES,DBYTES]=SEIZMOSIZE(DATA) also returns the expected
%    on-disk size of the header and data portions of the record(s) in DATA
%    in the arrays HBYTES and DBYTES.
%
%    Notes:
%
%    Examples:
%     % Calculate expected filesizes for SAC files in current directory:
%     seizmosize(readheader('*.SAC'))
%
%    See also: READHEADER, READDATA

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
%        Oct. 15, 2009 - speed boost
%        Jan. 30, 2010 - drop CHECKHEADER use, proper SEIZMO handling,
%                        versioninfo caching support, added extra outputs
%        Aug. 21, 2010 - nargchk fix
%        Feb. 11, 2011 - minor doc reformatting
%        Mar. 13, 2012 - use getheader improvements
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2012 at 23:55 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin));

% headers setup (check struct too)
[h,vi]=versioninfo(data);
h=[h.data];

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);
oldversioninfocache=versioninfo_cache(true);

% attempt size estimation
try
    % header info
    [npts,ncmp,iftype,leven]=getheader(data,...
        'npts','ncmp','iftype id','leven lgc');
    
    % filetype
    count=strcmpi(iftype,'itime')...
        +strcmpi(iftype,'ixy')+strcmpi(iftype,'ixyz')...
        +2*strcmpi(iftype,'irlim')+2*strcmpi(iftype,'iamph');

    % multi-component
    count=ncmp.*count;

    % uneven sampling
    count=count+strcmpi(leven,'false');

    % final tally
    hbytes=[h(vi).startbyte].';
    dbytes=count.*npts.*[h(vi).bytesize].';
    bytes=hbytes+dbytes;

    % toggle struct checking back
    seizmocheck_state(oldseizmocheckstate);
    versioninfo_cache(oldversioninfocache);
catch
    % toggle struct checking back
    seizmocheck_state(oldseizmocheckstate);
    versioninfo_cache(oldversioninfocache);
    
    % rethrow error
    error(lasterror);
end

end
