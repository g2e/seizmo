function [data]=removetrend(data)
%REMOVETREND    Remove linear fit from SEIZMO records
%
%    Usage:    data=removetrend(data)
%
%    Description:
%     DATA=REMOVETREND(DATA) removes the linear fit from SEIZMO records by
%     subtracting the best straight line fit to the data as determined by a
%     least squares inversion.  For multi-component records, each component
%     is dealt with separately.  It is highly recommended to combine this
%     command with any filtering operations to reduce edge effects that may
%     lead to poor data quality.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % 4th order lowpass butter filter with a passband corner of 10s
%     data=iirfilter(removetrend(data),'low','butter',1/10,4)
%
%    See also: REMOVEMEAN, REMOVEPOLYNOMIAL, GETPOLYNOMIAL, TAPER,
%              REMOVEDEADRECORDS, REMOVESPLINE, GETSPLINE, POLYFIT,
%              POLYVAL, DETREND

%     Version History:
%        Oct. 31, 2007 - initial version
%        Nov.  7, 2007 - doc update
%        Nov. 27, 2007 - RDRIFT added (removes both mean and trend)
%        Jan. 14, 2008 - handle uneven files
%        Feb. 12, 2008 - RTREND renamed to RSLOPE
%        Feb. 29, 2008 - SEISCHK support, handle uneven files better,
%                        RTREND removed (RSLOPE still around)
%        Mar.  4, 2008 - doc update, minor code cleaning
%        May  12, 2008 - fix dep* formula
%        June 12, 2008 - doc update, history added, renamed from RDRIFT to
%                        RTREND to match SAC, dropped RSLOPE
%        Oct.  3, 2008 - .dep & .ind
%        Nov. 22, 2008 - doc update, history fix, renamed from RTREND to
%                        REMOVETREND, one CHANGEHEADER call, better checks
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        June 24, 2009 - minor doc update
%        Dec.  6, 2009 - workaround: remove mean first for stability
%        Jan. 30, 2010 - seizmoverbose support, proper SEIZMO handling
%        Feb.  2, 2010 - versioninfo caching
%        Feb. 11, 2011 - mass nargchk fix, dropped versioninfo caching
%        Mar. 13, 2012 - doc update, seizmocheck fix, use getheader
%                        improvements
%        June  3, 2012 - no checkheader call, skip doubles conversion
%        Jan. 21, 2014 - minor doc fix
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

% attempt trend removal
try
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % header info
    leven=~strcmpi(getheader(data,'leven lgc'),'false');
    
    % detail message
    if(verbose)
        disp('Removing Trend from Record(s)');
        print_time_left(0,nrecs);
    end

    % remove trend and update header
    [depmen,depmin,depmax]=deal(nan(nrecs,1));
    for i=1:numel(data)
        % skip dataless
        if(isempty(data(i).dep))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

        % unevenly spaced
        if(~leven(i))
            for j=1:size(data(i).dep,2)
                data(i).dep(:,j)=data(i).dep(:,j) ...
                    -polyval(polyfit(double(data(i).ind),...
                    data(i).dep(:,j),1)-mean(data(i).dep(:,j)),...
                    double(data(i).ind));
            end
        else % evenly spaced
            for j=1:size(data(i).dep,2)
                data(i).dep(:,j)=detrend(...
                    data(i).dep(:,j)-mean(data(i).dep(:,j)));
            end
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
