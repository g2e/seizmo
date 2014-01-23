function [data]=removemean(data)
%REMOVEMEAN    Remove mean from SEIZMO records
%
%    Usage:    data=removemean(data)
%
%    Description:
%     DATA=REMOVEMEAN(DATA) removes the mean from SEIZMO records.  In the
%     case of multi-component records, each component has the mean removed.
%
%    Notes:
%     - REMOVEMEAN is useful for avoiding edge-effects in spectral
%       operations but REMOVETREND is probably a better option in this
%       case.
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % It is generally a good idea to remove the mean from records before
%     % performing any filtering operations to avoid edge effects:
%     plot1(squish(data,5))             % more ringing
%     plot1(squish(removemean(data),5)) % less ringing
%
%    See also: REMOVETREND, REMOVEPOLYNOMIAL, GETPOLYNOMIAL, TAPER,
%              REMOVEDEADRECORDS, REMOVESPLINE, GETSPLINE

%     Version History:
%        Oct. 31, 2007 - initial version
%        Nov.  7, 2007 - doc update
%        Nov. 27, 2007 - minor doc update
%        Feb. 29, 2008 - SEISCHK support
%        Mar.  4, 2008 - minor doc update
%        May  12, 2008 - fix dep* formula
%        June 12, 2008 - doc update, history added
%        Oct.  3, 2008 - .dep & .ind
%        Nov. 22, 2008 - doc update, rename from RMEAN to REMOVEMEAN
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        June 24, 2009 - minor doc fix
%        Jan. 29, 2010 - seizmoverbose support, proper SEIZMO handling
%        Feb. 11, 2011 - mass nargchk fix, mass seizmocheck fix
%        Apr.  3, 2012 - minor doc update
%        June  3, 2012 - skip doubles conversion
%        Jan. 21, 2014 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 21, 2014 at 15:05 GMT

% todo:

% check input
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt mean removal
try
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % detail message
    if(verbose)
        disp('Removing Mean from Record(s)');
        print_time_left(0,nrecs);
    end

    % remove mean and update header
    [depmen,depmin,depmax]=deal(nan(nrecs,1));
    for i=1:nrecs
        % skip dataless
        if(isempty(data(i).dep))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

        % loop through components
        for j=1:size(data(i).dep,2)
            data(i).dep(:,j)=data(i).dep(:,j)-mean(data(i).dep(:,j));
        end

        % adjust header
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % adjust header
    data=changeheader(data,...
        'depmen',depmen,'depmin',depmin,'depmax',depmax);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
