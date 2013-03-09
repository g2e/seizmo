function [data]=omegahilbert(data)
%OMEGAHILBERT    Hilbert Transform SEIZMO records in the frequency domain
%
%    Usage:    data=omegahilbert(data)
%
%    Description:
%     DATA=OMEGAHILBERT(DATA) returns the Hilbert transform of SEIZMO
%     spectral records in DATA.  The Hilbert transform can be thought of as
%     a -pi/2 phase shift.  This is useful for operations dealing with the
%     analytic signal of a record.
%
%    Notes:
%     - Because DFT zero-pads records this operation does not give the same
%       result of HILBRT (SEIZMO function) or HILBERT.
%
%    Examples:
%     % Compare H(u)(t) with u(t):
%     plot2([data(1) idft(omegahilbert(dft(data(1))))])
%
%    See also: HILBRT, ENVELOPE, OMEGASHIFT, OMEGADIVIDE, OMEGAMULTIPLY,
%              OMEGAANALYTIC, DFT, IDFT, OMEGAGAUSSIAN

%     Version History:
%        Feb.  4, 2012 - initial version
%        May  31, 2012 - minor doc update
%        Feb. 14, 2013 - use strcmpi for consistency
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 14, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt frequency domain hilbert transform
try
    % check headers
    data=checkheader(data,...
        'NONSPECTRAL_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR');

    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % retreive header info
    [npts,iftype]=getheader(data,'npts','iftype id');
    npts2=npts/2;
    
    % detail message
    if(verbose)
        disp('Hilbert Transforming Record(s) in the Frequency Domain');
        print_time_left(0,nrecs);
    end

    % loop through records
    [depmin,depmen,depmax]=deal(nan(nrecs,1));
    for i=1:nrecs
        % skip dataless
        if(isempty(data(i).dep))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);

        % hilbert transform
        if(strcmpi(iftype(i),'irlim'))
            % -i w>0, i w<0, 0 w=0
            data(i).dep(2:npts2(i),[1:2:end 2:2:end])=...
                [data(i).dep(2:npts2(i),2:2:end) ...
                -data(i).dep(2:npts2(i),1:2:end)];
            data(i).dep(npts2(i)+2:end,[1:2:end 2:2:end])=...
                [-data(i).dep(npts2(i)+2:end,2:2:end) ...
                +data(i).dep(npts2(i)+2:end,1:2:end)];
            data(i).dep([1 npts2(i)+1],:)=0;
        else % iamph
            data(i).dep(2:npts2(i),2:2:end)=...
                data(i).dep(2:npts2(i),2:2:end)-pi/2;
            data(i).dep(npts2(i)+2:end,2:2:end)=...
                data(i).dep(npts2(i)+2:end,2:2:end)+pi/2;
            data(i).dep([1 npts2(i)+1],:)=0;
        end

        % change class back
        data(i).dep=oclass(data(i).dep);

        % dep*
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % update header
    data=changeheader(data,...
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
