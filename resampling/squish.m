function [data]=squish(data,factor)
%SQUISH    Downsample SEIZMO records by an integer factor
%
%    Usage:    data=squish(data,factors)
%
%    Description:
%     SQUISH(DATA,FACTOR) downsamples SEIZMO records by an integer FACTOR
%     after implementing an anti-aliasing filter.  To avoid adding
%     significant numerical noise to the data, keep the decimation factor
%     below 13.  If a decimation factor greater than this is needed,
%     consider using a cascade of decimation factors (by giving an array of
%     factors for FACTOR).  For a usage case, see the examples below.
%     Strong amplitudes at or near the start and end of records will
%     introduce edge effects that can be avoided by first detrending and 
%     then tapering records.  Uses the Matlab function decimate (Signal 
%     Processing Toolbox).
%
%    Notes:
%
%    Header changes: DELTA, NPTS, E, DEPMEN, DEPMIN, DEPMAX
%
%    Examples: 
%     % To halve the samplerates of records in data:
%     data=squish(data,2);
%
%     % To cascade records to a samplerate 40 times lower:
%     data=squish(data,[5 8]);
%     % or:
%     data=squish(data,factor(40));
%
%    See also: STRETCH, SYNCRATES, INTERPOLATE

%     Version History:
%        Oct. 30, 2007 - initial version
%        Nov.  7, 2007 - doc update
%        Jan. 30, 2008 - better input checking and doc update
%        Feb. 23, 2008 - fixed bugs on multi-component files, improved
%                        input checks, added examples
%        Feb. 28, 2008 - seischk support
%        May  12, 2008 - dep* fix
%        June 15, 2008 - doc update
%        June 30, 2008 - single ch call, dataless support, handle 
%                        decimation to single point, doc update,
%                        .dep rather than .x, filetype checking
%        Nov. 23, 2008 - update for new name schema (now SQUISH)
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
    [b,delta,npts,ncmp]=getheader(data,'b','delta','npts','ncmp');
    
    % check factor
    if(~isnumeric(factor) || any(factor<1) || any(fix(factor)~=factor))
        error('seizmo:squish:badInput',...
            'FACTOR must be one or more positive integers!')
    end

    % number of factors
    factor=factor(:);
    nf=numel(factor);

    % overall factor
    of=prod(factor);
    delta=delta*of;
    
    % detail message
    if(verbose)
        disp('Decimating Record(s)');
        print_time_left(0,nrecs);
    end

    % decimate and update header
    e=nan(nrecs,1); depmen=e; depmin=e; depmax=e;
    for i=1:nrecs
        % dataless support
        if(~npts(i))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);

        % loop over components (because decimate can't handle arrays)
        for j=1:ncmp(i)
            temp=data(i).dep(:,j);
            % loop over factors
            for k=1:nf
                temp=decimate(temp,factor(k));
            end
            % preallocate after we know length
            if(j==1)
                npts(i)=size(temp,1);
                e(i)=b(i)+(npts(i)-1)*delta(i);
                save=zeros(npts(i),ncmp(i));
            end
            save(:,j)=temp;
        end

        % change class back
        data(i).dep=oclass(save);

        % dep* stats
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % update header
    data=changeheader(data,'npts',npts,'delta',delta,'e',e,...
        'depmin',depmin,'depmen',depmen,'depmax',depmax);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
