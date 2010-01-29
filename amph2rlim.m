function [data]=amph2rlim(data)
%AMPH2RLIM    Convert SEIZMO spectral records from AMPH to RLIM
%
%    Usage:    data=amph2rlim(data)
%
%    Description: AMPH2RLIM(DATA) converts SEIZMO amplitude-phase records 
%     to real-imaginary records.  This is particularly useful when
%     performing basic operations on spectral records which would otherwise
%     require treating the amplitude and phase components separately.
%     Records in DATA must be of the spectral variety.  Real-imaginary
%     records are not altered.
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
    iamph=find(amph); namph=numel(iamph);

    % records must be spectral
    if(any(~amph & ~rlim))
        error('seizmo:amph2rlim:illegalOperation',...
            ['Record(s):\n' sprintf('%d ',find(~amph & ~rlim)) ...
            '\nIllegal operation on non-spectral record(s)!']);
    end
    
    % detail message
    if(verbose)
        disp('Converting AMPH Record(s) to RLIM');
        print_time_left(0,nrecs);
    end
    
    % quick exit if all rlim
    if(namph==0)
        % detail message
        if(verbose)
            print_time_left(nrecs,nrecs);
        end
        return;
    end
    
    % loop through records
    depmen=nan(namph,1); depmin=depmen; depmax=depmen;
    for i=1:namph
        k=iamph(i);

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
        data(i).dep=double(data(k).dep);
        temp=data(k).dep(:,1:2:end).*exp(j*data(k).dep(:,2:2:end));
        data(k).dep(:,1:2:end)=real(temp);
        data(k).dep(:,2:2:end)=imag(temp);
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
    data(amph)=changeheader(data(amph),'iftype','irlim',...
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
