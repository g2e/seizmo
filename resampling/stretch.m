function [data]=stretch(data,factor)
%STRETCH    Upsample SEIZMO records by an integer factor
%
%    Usage:    data=stretch(data,factor)
%
%    Description:
%     STRETCH(DATA,FACTOR) upsamples SEIZMO records by an integer FACTOR.
%     Anti-aliasing is not an issue so there is no limit imposed on the
%     stretch factor.  Cascades are allowed regardless.  Uses the Matlab
%     function interp (Signal Processing Toolbox).  Detrending and tapering
%     records beforehand will help to suppress edge effects.
%
%    Notes:
%
%    Header changes: DELTA, NPTS, DEPMEN, DEPMIN, DEPMAX
%
%    Examples: 
%     % To double samplerates:
%     data=stretch(data,2);
%
%     % To cascade to a samplerate 40 times higher:
%     data=stretch(data,[5 8]);
%
%    See also: SQUISH, SYNCRATES, INTERPOLATE

%     Version History:
%        Oct. 31, 2007 - initial version
%        Nov.  7, 2007 - doc update
%        Jan. 30, 2008 - doc update, code cleaning
%        Mar.  4, 2008 - major code cleaning, fix extrapolation
%        May  12, 2008 - fix dep* formula
%        June 15, 2008 - doc update, history added
%        Nov. 23, 2008 - doc update, history fixed, better checking, single
%                        changeheader call, update for new name schema,
%                        .dep rather than .x
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Jan. 30, 2010 - proper SEIZMO handling, seizmoverbose support,
%                        improved error messages
%        Feb.  2, 2010 - versioninfo caching
%        Feb. 11, 2011 - mass nargchk fix, dropped versioninfo caching
%        Mar. 13, 2012 - doc update, seizmocheck fix, better checkheader
%                        usage
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2012 at 15:05 GMT

% todo:

% check inputs
error(nargchk(2,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% empty factor
if(isempty(factor)); return; end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt convolution
try
    % check headers
    data=checkheader(data,...
        'FALSE_LEVEN','ERROR',...
        'NONTIME_IFTYPE','ERROR');

    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % get header info
    [delta,npts,ncmp]=getheader(data,'delta','npts','ncmp');

    % check factor
    if(~isnumeric(factor) || any(factor<1) || any(fix(factor)~=factor))
        error('seizmo:stretch:badInput',...
            'FACTOR must be one or more positive integers!');
    end

    % number of factors
    factor=factor(:);
    nf=numel(factor);

    % overall factor
    of=prod(factor);
    delta=delta*of;

    % detail message
    if(verbose)
        disp('Upsampling Record(s)');
        print_time_left(0,nrecs);
    end

    % decimate and update header
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
    for i=1:nrecs
        % skip dataless
        if(~npts(i))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).x=double(data(i).dep);

        % loop over components (b/c interp can't handle arrays)
        npts(i)=(npts(i)-1)*of+1;
        save=zeros(npts(i),ncmp(i));
        for j=1:ncmp(i)
            temp=data(i).dep(:,j);
            % loop over factors
            for k=1:nf
                % stretch (filter length 4, cutoff freq is nyquist)
                temp=interp(temp,factor(k),4,1);
            end

            % Matlab extrapolates past the end, so here
            % we truncate the extrapolated values off
            save(:,j)=temp(1:npts(i),1);
        end

        % change class back
        data(i).dep=oclass(save);

        % get dep*
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));

        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % update header
    data=changeheader(data,'npts',npts,'delta',delta,...
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
