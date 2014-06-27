function [data]=amph2rlim(data)
%AMPH2RLIM    Convert SEIZMO spectral records from AMPH to RLIM
%
%    Usage:    data=amph2rlim(data)
%
%    Description:
%     AMPH2RLIM(DATA) converts SEIZMO amplitude-phase records to
%     real-imaginary records.  This is particularly useful when performing
%     basic operations on spectral records which would otherwise require
%     treating the amplitude and phase components separately.
%     Real-imaginary and other record filetypes are not altered.
%
%    Notes:
%
%    Header changes: IFTYPE, DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % It is easier to multiply records by a constant in the frequency
%     % domain if they are in real-imaginary format:
%     data=amph2rlim(data);
%     data=multiply(data,3);
%     data=rlim2amph(data);
%
%    See also: RLIM2AMPH, DFT, IDFT

%     Version History:
%        June 11, 2008 - initial version
%        June 20, 2008 - minor doc update
%        June 28, 2008 - fixed call to ch, removed option,
%                        doc update, .dep rather than .x
%        July 19, 2008 - dataless support, updates DEP* fields
%        Oct.  7, 2008 - minor code cleaning
%        Nov. 22, 2008 - update for new name schema
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Oct. 21, 2009 - only touches amph (maybe a bit faster)
%        Dec.  4, 2009 - fixed IFTYPE bug, handle no amph case
%        Jan. 26, 2010 - seizmoverbose support
%        Feb.  2, 2010 - versioninfo caching (required some code changes)
%        Mar.  8, 2010 - versioninfo caching dropped
%        Aug. 16, 2010 - nargchk fix, better checkheader utilization
%        Dec. 21, 2011 - doc update, changed example (it was bad)
%        Mar. 13, 2012 - use getheader improvements
%        June 12, 2014 - remove checkheader call, allow non-spectral file
%                        passthrough, no verbose messages
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2014 at 12:45 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% retreive header info
try
    iftype=getheader(data,'iftype id');
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% number of records
nrecs=numel(data);

% find spectral
amph=strcmpi(iftype,'iamph');

% anything amph? if no, skip
if(~any(amph)); return; end

% loop through records
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % only bother with amph
    if(amph(i))
        % convert amph
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);
        temp=data(i).dep(:,1:2:end).*exp(1j*data(i).dep(:,2:2:end));
        data(i).dep(:,1:2:end)=real(temp);
        data(i).dep(:,2:2:end)=imag(temp);
        data(i).dep=oclass(data(i).dep);
        
        % dep*
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
    end
end

% update filetype
data(amph)=changeheader(data(amph),'iftype','irlim',...
    'depmax',depmax(amph),'depmin',depmin(amph),'depmen',depmen(amph));

end
