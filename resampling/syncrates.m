function [data]=syncrates(data,sr,tol)
%SYNCRATES    Resample SEIZMO records to a common sample rate
%
%    Usage:    data=syncrates(data,sr)
%              data=syncrates(data,sr,tol)
%
%    Description:
%     SYNCRATES(DATA,SR) syncronizes the sample rates of SEIZMO records in
%     DATA to the sample rate SR.  A fir filter is used to avoid aliasing
%     issues, but this can cause edge effects if the records deviate from
%     zero strongly at the start/end of the record.  Typically using
%     REMOVETREND and TAPER on records beforehand helps to limit the edge
%     effects.  Uses the Matlab function RESAMPLE (Signal Processing
%     Toolbox).
%
%     SYNCRATES(DATA,SR,TOL) specifies the maximum tolerance TOL that the
%     fraction of 2 small integers must match the ratio of the old and new
%     sample rates of a record.  The integers specify the upsampling and
%     downsampling portions of the resampling operation.  See function RRAT
%     for more details.  The default TOL is 1e-6:
%       abs(Up/Down - New/Old) / (New/Old) <= 1e-6
%
%    Notes:
%     - requires evenly sampled data (use INTERPOLATE for uneven data)
%
%    Header Changes: DELTA, NPTS, DEPMEN, DEPMIN, DEPMAX, E
%
%    Examples:
%     % Change all records to 5 samples per second:
%     data=syncrates(data,5);
%
%    See also: INTERPOLATE, IIRFILTER, SQUISH, STRETCH

%     Version History:
%        Feb. 16, 2008 - initial version
%        Feb. 23, 2008 - minor fix
%        Feb. 26, 2008 - minor fix
%        Mar.  4, 2008 - doc update, better checks
%        June 16, 2008 - doc update, rejects spectral records
%                        fixed no header update bug, and changed the way
%                        records are resampled when Matlab's resample fails
%        Nov. 22, 2008 - update for new name schema (now SYNCRATES), allow
%                        resampling spectral records, disallow resampling
%                        xyz records
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        June 25, 2009 - update for RECORD2MAT/MAT2RECORD, process records
%                        individually, minor bug fix for rare case when
%                        resample does not work
%        Nov. 26, 2009 - document RAT bug, alter RAT call slightly to force
%                        better accuracy of the resampling operation, add
%                        TOL argument, fix NPTS handling
%        Jan. 30, 2010 - proper SEIZMO handling, seizmoverbose support,
%                        improved error messages, use RRAT (fixed RAT)
%        Feb.  3, 2010 - make sure checkheader is skipped for trouble case
%        Feb. 11, 2011 - mass nargchk fix, mass seizmocheck fix
%        Mar. 13, 2012 - disallow spectral syncing, doc update, better
%                        checkheader usage
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    data=checkheader(data,...
        'FALSE_LEVEN','ERROR',...
        'NONTIME_IFTYPE','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% attept synchronization of sample rates
try
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);

    % get header info
    [delta,b]=getheader(data,'delta','b');

    % check rate
    if(~isnumeric(sr) || ~isscalar(sr) || sr<=0)
        error('seizmo:syncrates:badInput',...
            'SR must be a positive numeric scalar!');
    end

    % check tol
    if(nargin==2 || isempty(tol)); tol=1e-6; end
    if(~isscalar(tol) || ~isreal(tol))
        error('seizmo:syncrates:badInput',...
            'TOL must be a scalar real!');
    end

    % find fraction numerator/denominator of sampling
    % rate ratio expressed as integers
    [n,d]=rrat(delta*sr,tol);
    
    % detail message
    if(verbose)
        disp('Syncing Sample Rates of Record(s)');
        print_time_left(0,nrecs);
    end

    % loop over every record
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen; npts=depmen;
    for i=1:nrecs
        % skip dataless
        if(isempty(data(i).dep))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

        % save class
        oclass=str2func(class(data(i).dep));

        % try resample
        try
            data(i).dep=oclass(resample(double(data(i).dep),n(i),d(i)));
        catch
            disp(['Had some trouble resampling record: ' num2str(i)])
            % interpolate to a sample rate at least 4 times the current
            % AND the desired sample rates to preserve the data and allow
            % for decimation.  here we find the lowest multiple of the
            % desired rate that is 4x the current rate...hopefully that
            % integer is fairly factorable...if not, the function crashes!
            lm=max([4 ceil(4/(delta(i)*sr))]);

            % interpolate temporarily to higher rate
            data(i)=interpolate(data(i),lm*sr,'spline');

            % now resample
            data(i).dep=oclass(resample(double(data(i).dep),1,lm));
        end

        % get npts
        npts(i)=size(data(i).dep,1);

        % get dep*
        if(npts(i))
            depmen(i)=nanmean(data(i).dep(:));
            depmin(i)=min(data(i).dep(:));
            depmax(i)=max(data(i).dep(:));
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % update header
    data=changeheader(data,'delta',1/sr,'npts',npts,'e',b+(npts-1)./sr,...
        'depmen',depmen,'depmin',depmin,'depmax',depmax);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror);
end

end
