function [xp,xi]=getpeaks(x,varargin)
%GETPEAKS    Returns peak info of a data vector
%
%    Usage:    [xp,xi]=getpeaks(x)
%              [xp,xi]=getpeaks(...,'npeaks',n,...)
%              [xp,xi]=getpeaks(...,'spacing',s,...)
%              [xp,xi]=getpeaks(...,'adjacent',a,...)
%              [xp,xi]=getpeaks(...,'highest',h,...)
%              [xp,xi]=getpeaks(...,'fast',f,...)
%
%    Description:
%     [XP,XI]=GETPEAKS(X) returns the highest peak in X.  X must be a
%     numeric vector.  XP is the value of the peak and XI is the index of
%     the peak in X.  This is the same as [XP,XI]=MAX(X).
%
%     [XP,XI]=GETPEAKS(...,'NPEAKS',N,...) returns the highest N peaks
%     where peaks are required to have a positive difference from both of
%     their neighboring points.  Peaks are sorted by the value of the peak
%     not by the height from their neighbors.
%
%     [XP,XI]=GETPEAKS(...,'SPACING',S,...) requires peaks to be separated
%     by at least S samples.  The default spacing is S=1 but all peaks will
%     be separated by at least 1 sample anyways due to the nature of the
%     peak finding.
%
%     [XP,XI]=GETPEAKS(...,'ADJACENT',A,...) also returns the points
%     adjacent to the selected peaks within A samples.  The default is A=0.
%     Setting A=1 would grab the points directly adjacent to the peaks as
%     well making XP & XI an Nx3 array.  Values and indices go from -A to A
%     for each peak (row).
%
%     [XP,XI]=GETPEAKS(...,'HIGHEST',H,...) returns the highest peaks if
%     H=TRUE or the lowest peaks if H=FALSE.  The default is H=TRUE.
%
%     [XP,XI]=GETPEAKS(...,'FAST',F,...) uses the fast peak algorithm if
%     F=TRUE and N=1 (the defaults).  This is substantially faster because
%     it just calls MAX.  Setting N=2+ will ignore this parameter.
%
%    Notes:
%
%    Examples:
%     % Get the 5 deepest troughs of some random data:
%     x=rand(100,1);
%     [xp,xi]=getpeaks(-x,'n',5);
%     figure; plot(x); hold on; plot(xi,-xp,'ro'); hold off;
%     
%     % Highlight the 4 highest peaks and their neighboring points:
%     x=rand(100,1);
%     [xp,xi]=getpeaks(x,'n',4,'a',1);
%     figure; plot(x); hold on;
%     plot(xi',xp','r'); plot(xi(:,2),xp(:,2),'ro'); hold off;
%
%    See also: CORRELATE, FINDPEAKS

%     Version History:
%        Oct. 18, 2012 - initial version
%        Oct. 19, 2012 - tweaks for speed
%        Aug. 16, 2013 - bugfixes for no peak detection
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 16, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check x
if(~isvector(x) || ~isnumeric(x))
    error('seizmo:getpeaks:badInput',...
        'X must be a numeric vector!');
end

% default/check/get options
n=1; s=1; a=0; h=true; f=true;
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:getpeaks:badInput',...
        'Parameters must be specified with strings!');
end
valid=['npeaks  '; 'spacing '; 'adjacent'; 'highest ';'fast    '];
for i=1:2:nargin-1
    switch strmatch(lower(varargin{i}),valid)
        case 1 % npeaks
            n=varargin{i+1};
            if(~isscalar(n) || ~isnumeric(n) || ~isreal(n) || n<=0 ...
                    || n~=fix(n))
                error('seizmo:getpeaks:badInput',...
                    'NPEAKS must be an integer >=1 !');
            end
        case 2 % spacing
            s=varargin{i+1};
            if(~isscalar(s) || ~isnumeric(s) || ~isreal(s) || s<=0 ...
                    || s~=fix(s))
                error('seizmo:getpeaks:badInput',...
                    'SPACING must be an integer >=1 !');
            end
        case 3 % adjacent
            a=varargin{i+1};
            if(~isscalar(a) || ~isnumeric(a) || ~isreal(a) || a<=0 ...
                    || a~=fix(a))
                error('seizmo:getpeaks:badInput',...
                    'ADJACENT must be an integer >=0 !');
            end
        case 4 % highest
            h=varargin{i+1};
            if(~isscalar(h) || ~islogical(h))
                error('seizmo:getpeaks:badInput',...
                    'HIGHEST must be TRUE or FALSE!');
            end
        case 5 % fast
            f=varargin{i+1};
            if(~isscalar(f) || ~islogical(f))
                error('seizmo:getpeaks:badInput',...
                    'FAST must be TRUE or FALSE!');
            end
        otherwise
            error('seizmo:getpeaks:badInput',...
                'Unknown Parameter: %s !',varargin{i});
    end
end

% fast?
if(f && n==1)
    [xp,xi]=max(x);
else
    % force column vector
    x=x(:);
    
    % all peaks toward positive infinity
    xi=find(diff([-inf; x])>0 & diff([x; -inf])<0);
    xp=x(xi);
    
    % sort by peak height
    if(h); [xp,idx]=sort(xp,'descend');
    else [xp,idx]=sort(xp,'ascend');
    end
    xi=xi(idx);
    
    % eliminate by spacing
    if(s>1)
        peak=1;
        while(peak<min(n+1,numel(xi)))
            % eliminate those within s units of current peak
            xi([false(peak,1); abs(xi(peak)-xi(peak+1:end))<=s])=[];
            peak=peak+1;
        end
    end
    
    % pad/clip if necessary
    xi(max(end,1):n)=nan;
    xi(n+1:end)=[];
    xp(max(end,1):n)=nan;
    xp(n+1:end)=[];
end

% grab adjacent points
if(a>0)
    adj=-a:a;
    xi=xi(:,ones(1,2*a+1))+adj(ones(n,1),:);
    xp=nan(size(xi));
    xp(xi>0 & xi<=numel(x))=x(xi(xi>0 & xi<=numel(x)));
end

end
