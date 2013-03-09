function [data]=omegagaussian(data,fc,a)
%OMEGAGAUSSIAN    Gaussian filter of SEIZMO records in the frequency domain
%
%    Usage:    data=omegagaussian(data,fc)
%              data=omegagaussian(data,fc,a)
%
%    Description:
%     DATA=OMEGAGAUSSIAN(DATA,FC) applies a gaussian filter to the SEIZMO
%     spectral records in DATA.  The gaussian filter is centered at FC (in
%     Hertz) and is defined as:
%
%                                2
%                     / (f-fc) \
%                 -A |  ______  |
%                     \   fc   /
%     (1) G(f) = e
%
%     (2) D(f) = D(f)G(f)
%
%     where G(f) is the gaussian filter, D(f) is the data, A is a tunable
%     parameter controlling the frequency vs time resolution, f are the
%     frequencies, & fc is the center frequency of the gaussian.  A by
%     default is 100.
%
%     DATA=OMEGAGAUSSIAN(DATA,FC,A) adjusts the resolution parameter A
%     defined above.  Higher A gives better frequency resolution while
%     lower A gives better time resolution.  A is related to the width of
%     the gaussian filter at a specified dropoff from the peak by:
%
%                B
%      (3) A = _____
%                 2
%                W
%
%      where B is the dropoff fraction expressed as the natural logarithm
%      of the height at the center frequency to the dropoff height and W is
%      the halfwidth of the gaussian at the desired dropoff height
%      expressed as a fraction of FC.  See the Notes & Examples sections
%      for further explanation.
%
%    Notes:
%     - The following paper describes the algorithms used here:
%        Dziewonski, Bloch, & Landisman 1969, BSSA Vol 59, No 1, pp 427-444
%
%    Examples:
%     % Gaussian filter centered at .01Hz with a width of .005Hz when the
%     % gaussian has dropped off to 1/100th of what it was at the center
%     % frequency (it is 1 at the center frequency):
%     data=omegagaussian(data,.01,log(100)/(.005/2/.01)^2);
%
%    See also: OMEGAHILBERT, OMEGASHIFT, OMEGADIVIDE, OMEGAMULTIPLY,
%              OMEGAANALYTIC, DFT, IDFT, IIRFILTER, TAPERFUN, FTAN

%     Version History:
%        Feb.  5, 2012 - initial version
%        Aug.  2, 2012 - minor doc update
%        Feb. 14, 2013 - use strcmpi for consistency
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 14, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt frequency domain gaussian filter
try
    % check headers
    data=checkheader(data,...
        'NONSPECTRAL_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR');

    % verbosity
    verbose=seizmoverbose;
    
    % check fc & a
    if(nargin<3 || isempty(a)); a=100; end
    if(~isnumeric(fc) || ~isreal(fc) || any(fc<0))
        error('seizmo:omegagaussian:badInput',...
            'FC must be a positive real value in Hz!');
    elseif(~isnumeric(a) || ~isreal(a) || any(a<0))
        error('seizmo:omegagaussian:badInput',...
            'A must be a positive real value!');
    end
    [data,fc,a]=expandscalars(data,fc,a);
    
    % number of records
    nrecs=numel(data);
    
    % retreive header info
    [delta,npts,ncmp,iftype]=getheader(data,...
        'delta','npts','ncmp','iftype id');
    npts2=npts/2;
    
    % detail message
    if(verbose)
        disp('Gaussian Filtering Record(s) in the Frequency Domain');
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

        % gaussian filter
        w=delta(i)*[0:1:npts2(i) npts2(i)-1:-1:1].';
        g=exp(-a(i)*((w-fc(i))/fc(i)).^2);
        if(strcmpi(iftype(i),'irlim'))
            data(i).dep=data(i).dep.*g(:,ones(1,ncmp(i)));
        else % iamph
            data(i).dep(:,1:2:end)=data(i).dep(:,1:2:end)...
                .*g(:,ones(1,ncmp(i)));
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
