function [data]=rlim2amph(data)
%RLIM2AMPH    Convert SEIZMO spectral records from RLIM to AMPH
%
%    Usage:    data=rlim2amph(data)
%
%    Description: RLIM2AMPH(DATA) converts SEIZMO real-imaginary records 
%     to amplitude-phase records.  This is particularly useful for 
%     switching between the formats when performing basic operations on 
%     spectral records would require separating the amplitude and phase
%     components.  Records in DATA must be of the spectral variety.  
%     Amplitude-phase records are not altered.
%
%    Notes:
%
%    Header changes: IFTYPE, DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     To simply multiply two records in the frequency domain, they must be
%     converted to real-imaginary first:
%      data=amph2rlim(data)
%      data=multiplyrecords(data(1),data(2))
%      data=rlim2amph(data)
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2010 at 18:50 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt conversion
try
    % check header
    data=checkheader(data);
    
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);
    
    % retreive header info
    iftype=getenumid(data,'iftype');

    % find spectral
    amph=strcmpi(iftype,'iamph');
    rlim=strcmpi(iftype,'irlim');
    irlim=find(rlim); nrlim=numel(irlim);

    % records must be spectral
    if(any(~amph & ~rlim))
        error('seizmo:rlim2amph:illegalOperation',...
            ['Record(s):\n' sprintf('%d ',find(~amph & ~rlim)) ...
            '\nIllegal operation on non-spectral record(s)!']);
    end
    
    % detail message
    if(verbose)
        disp('Converting RLIM Record(s) to AMPH');
        print_time_left(0,nrecs);
    end
    
    % quick exit if all amph
    if(nrlim==0)
        % detail message
        if(verbose)
            print_time_left(nrecs,nrecs);
        end
        return;
    end

    % loop through records
    depmen=nan(nrlim,1); depmin=depmen; depmax=depmen;
    for i=1:nrlim
        k=irlim(i);

        % skip dataless
        if(isempty(data(k).dep))
            % detail message
            if(verbose)
                print_time_left(k,nrecs);
            end
            continue;
        end

        % convert
        oclass=str2func(class(data(k).dep));
        data(k).dep=double(data(k).dep);
        temp=complex(data(k).dep(:,1:2:end),data(k).dep(:,2:2:end));
        data(k).dep(:,1:2:end)=abs(temp);
        data(k).dep(:,2:2:end)=angle(temp);
        data(k).dep=oclass(data(k).dep);

        % dep*
        depmen(i)=mean(data(k).dep(:));
        depmin(i)=min(data(k).dep(:));
        depmax(i)=max(data(k).dep(:));
        
        % detail message
        if(verbose)
            print_time_left(k,nrecs);
        end
    end
    
    % detail message
    if(verbose && k~=nrecs)
        print_time_left(nrecs,nrecs);
    end

    % update filetype
    data(rlim)=changeheader(data(rlim),'iftype','iamph',...
        'depmax',depmax,'depmin',depmin,'depmen',depmen);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
