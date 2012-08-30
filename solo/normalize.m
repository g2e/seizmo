function [data,scale]=normalize(data)
%NORMALIZE    Normalizes SEIZMO records
%
%    Usage:    [data,scale]=normalize(data)
%
%    Description:
%     [DATA,SCALE]=NORMALIZE(DATA) scales the amplitudes of SEIZMO records
%     to the range of -1 to 1.  SCALE contains the normalization factors.
%     All components for a record are scaled by the same factor, which is
%     the maximum amplitude found assuming the components are orthogonal.
%     The maximum amplitude is just:
%                        _________________
%                       /
%       amp_max=max \  / cmp1^2+cmp2^2+...
%                    \/
%
%     Use header fields DEPMIN and DEPMAX to get the single largest data
%     value.
%
%    Notes:
%     - this formulation normalizes real-imaginary spectral records too
%
%    Header changes: DEPMIN, DEPMAX, DEPMEN
%
%    Examples:
%     % Display record overlay with records scaled:
%     plot2(normalize(data))
%
%    See also: MULTIPLY, GETNORM

%     Version History:
%        Feb. 12, 2008 - initial version
%        Feb. 23, 2008 - bug fix
%        Feb. 28, 2008 - real bug fix
%        Mar.  4, 2008 - minor doc update
%        Apr. 17, 2008 - minor doc update
%        June 16, 2008 - doc update
%        Nov. 23, 2008 - renamed from NRM to NORMALIZE, update for new name
%                        schema, .dep over .x, one changeheader call
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Jan. 30, 2010 - proper SEIZMO handling, seizmoverbose support
%        Feb. 11, 2011 - mass nargchk fix, mass seizmocheck fix
%        Apr. 14, 2011 - fix for record of all zeros returning nans
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt normalization
try
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % detail message
    if(verbose)
        disp('Normalizing Record(s)')
        print_time_left(0,nrecs);
    end

    % normalize data
    scale=nan(nrecs,1); depmen=scale; depmin=scale; depmax=scale;
    for i=1:nrecs
        % skip dataless
        if(isempty(data(i).dep))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

        % get class and change to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);

        % get norm
        scale(i)=max(sqrt(sum((data(i).dep).^2,2)));

        % scale data and change class back
        if(scale(i))
            data(i).dep=oclass(data(i).dep/scale(i));
        end

        % get dep*
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
    
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % update header
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
