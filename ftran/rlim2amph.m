function [data]=rlim2amph(data)
%RLIM2AMPH    Convert SEIZMO spectral records from RLIM to AMPH
%
%    Usage:    data=rlim2amph(data)
%
%    Description:
%     RLIM2AMPH(DATA) converts SEIZMO real-imaginary records to
%     amplitude-phase records.  This is particularly useful for switching
%     between the formats when performing basic operations on spectral
%     records would require separating the amplitude and phase components.
%     Records in DATA must be of the spectral variety.  Amplitude-phase
%     records are not altered.
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
%    See also: AMPH2RLIM, DFT, IDFT

%     Version History:
%        June 11, 2008 - initial version
%        July 19, 2008 - removed option, single call to ch, dataless
%                        support, updates DEP* fields, .dep rather than .x,
%                        doc update
%        Oct.  7, 2008 - minor code cleaning
%        Nov. 22, 2008 - update for new name schema
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Oct. 21, 2009 - only touches rlim (maybe a bit faster)
%        Dec.  4, 2009 - handle no rlim case
%        Jan. 26, 2010 - seizmoverbose support
%        Feb.  2, 2010 - versioninfo caching (required some code changes)
%        Mar.  8, 2010 - versioninfo caching dropped
%        Apr.  9, 2010 - minor bug fix
%        Feb. 11, 2011 - mass nargchk fix
%        Dec. 21, 2011 - doc update, changed example (it was bad)
%        Mar. 13, 2012 - use getheader improvements
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2012 at 12:45 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt conversion
try
    % check header
    data=checkheader(data,...
        'NONSPECTRAL_IFTYPE','ERROR');
    
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);
    
    % retreive header info
    iftype=getheader(data,'iftype id');

    % find spectral
    rlim=strcmpi(iftype,'irlim');
    
    % detail message
    if(verbose)
        disp('Converting RLIM Record(s) to AMPH');
        print_time_left(0,nrecs);
    end

    % loop through records
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
    for i=1:nrecs
        % skip dataless
        if(isempty(data(i).dep))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end
        
        % convert rlim
        if(rlim(i))
            oclass=str2func(class(data(i).dep));
            data(i).dep=double(data(i).dep);
            temp=complex(data(i).dep(:,1:2:end),data(i).dep(:,2:2:end));
            data(i).dep(:,1:2:end)=abs(temp);
            data(i).dep(:,2:2:end)=angle(temp);
            data(i).dep=oclass(data(i).dep);
        end

        % dep*
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % update filetype
    data=changeheader(data,'iftype','iamph',...
        'depmax',depmax,'depmin',depmin,'depmen',depmen);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
