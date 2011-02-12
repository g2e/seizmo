function [times,n]=getarrival(data,phase)
%GETARRIVAL    Returns stored phase arrival time from SEIZMO data header
%
%    Usage:    times=getarrival(data,phase)
%              [times,n]=getarrival(data,phase)
%
%    Description: GETARRIVAL(DATA,PHASE) searches 'kt(n)' header fields in
%     the SEIZMO structure DATA for the specified phase PHASE.  If found,
%     the matching 't(n)' value is returned.  If not, NaN is returned. In
%     case of multiple entries for the same phase, only the first match
%     found for each record is returned - lower index has preference.  Note
%     that the returned time is based on the reference time and is not
%     relative to the origin (see example below).  Phase may be a string
%     like 'Pn' or a cell array of strings like {'P' 'Pdiff'}.  Note that
%     in the case of a list of phases only the first match is returned for
%     each record.
%
%     [TIMES,N]=GETARRIVAL(DATA,PHASE) also returns an index array that
%     indicates the header field from which the time was found for each
%     record.  So [3 2 1] would mean record 1's time came from t3, record
%     2's time from t2 and record 3's from t1.  This is particularly useful
%     for setting iztype when combining getarrival with timeshift.
%     
%    Notes:
%     - NAME OF PHASE IS CASE SENSITIVE!
%     - Returned times are relative to the reference time!
%
%    Examples:
%     Ptimes=getarrival(data,'P');
%     sSKStimes=getarrival(data,'sSKS')
%
%     Get arrival time that is relative to origin time:
%      Ptimes=getarrival(data,'P')-getheader(data,'o');
%
%    See also: QUICKSNR, GETHEADER, ADDARRIVALS, TIMESHIFT

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 29, 2008 - uses SEISCHK now
%        Mar.  4, 2008 - minor doc update
%        Oct.  8, 2008 - doc update, add history
%        Nov. 16, 2008 - doc update, history fix, rename from PULLARR to
%                        GETARRIVAL, minor code clean
%        Nov. 24, 2008 - minor code cleaning
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        June 29, 2009 - doc update
%        June 30, 2009 - second output: t index
%        Dec.  8, 2009 - minor doc fix
%        Jan. 29, 2010 - dropped most checks, seizmoverbose support
%        Feb. 24, 2010 - 1 warning for unfound arrivals
%        Mar.  1, 2010 - allow multi-phase search
%        Mar. 16, 2010 - fixed error message
%        Feb. 11, 2011 - mass nargchk fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% verbosity
verbose=seizmoverbose;

% number of records
nrecs=numel(data);

% grab header values
[kt,t]=getheader(data,'kt','t');

% remove spaces from kt
kt=strtrim(kt);

% check phase
if(ischar(phase)); phase=cellstr(phase); end
if(~iscellstr(phase))
    error('seizmo:getarrival:badInput',...
        'PHASE must be a string or cellstr of phases!');
end
nphases=numel(phase);

% detail message
if(verbose)
    disp('Finding Phase Info in Record Header(s)');
    print_time_left(0,nrecs);
end

% do operations individually
times=nan(nrecs,1); n=times;
for i=1:nrecs
    % loop over each phase until something is found
    for j=1:nphases
        % find first match
        pos=find(strcmp(phase{j},kt(i,:)),1);
        
        % check for something
        if(~isempty(pos)); break; end
    end
    
    % check for failure
    if(isempty(pos))
        % detail message
        if(verbose); print_time_left(i,nrecs); end
        continue;
    end
    
    % add time
    times(i)=t(i,pos);
    n(i)=pos-1;
    
    % detail message
    if(verbose); print_time_left(i,nrecs); end
end

% throw warning for missed arrivals
if(any(isnan(n)))
    if(all(isnan(n)))
        tmp=['ALL OF THEM: 1 to ' num2str(numel(n))];
    elseif(sum(isnan(n))>20)
        tmp=[sprintf('%d ',find(isnan(n),10,'first')) ...
            '... ' sprintf('%d ',find(isnan(n),10,'last'))];
    else
        tmp=sprintf('%d ',find(isnan(n)));
    end
    warning('seizmo:getarrival:noPhase',...
        ['Could not find phase(s):\n ' sprintf('%s ',phase{:}) ...
        '\nRecord(s):\n ' tmp]);
end

end
