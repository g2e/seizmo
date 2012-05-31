function [pts]=taperfun(type,pts,limits,varargin)
%TAPERFUN    Returns a taper as specified
%
%    Usage:    tpr=taperfun(type,t,lim)
%              tpr=taperfun(type,t,lim,opt)
%
%    Description:
%     TPR=TAPERFUN(TYPE,T,LIM) returns a leading taper of type TYPE.  The
%     taper goes from 0 to 1 within the time limits defined in LIM as
%     [lower_bound upper_bound].  T is an array of time points for which
%     corresponding taper values are returned as TPR.  Points below the
%     lower bound return 0, while points above the upper bound give 1.  If
%     LIM is specified as [upper_bound lower_bound] then a trailing taper
%     is used with points below the lower bound returning 1 and points
%     above the upper bound give 0.
%
%     TPR=TAPERFUN(TYPE,T,LIM,OPT) passes taper option OPT.  Only
%     'chebwin', 'gausswin', 'kaiser', and 'tukeywin' have options.  See
%     their doc pages for more specific info.  Setting OPT to [] or nan
%     will use the default value.  Defaults:
%      chebwin   100
%      gausswin  2.5
%      kaiser    0.5
%      tukeywin  0.5
%
%    Notes:
%     - Valid Types:
%        BARTHANNWIN, BARTLETT, BLACKMAN, BLACKMANHARRIS, BOHMANWIN,
%        CHEBWIN, COSINE, FLATTOPWIN, GAUSSWIN, HAMMING, HANN, KAISER,
%        NUTTALLWIN, PARZENWIN, RECTWIN, TUKEYWIN
%
%    Examples:
%     % Plot random values for a gaussian taper:
%     t=sort(rand(100,1));
%     plot(t,taperfun('gausswin',t,[0 .5],4),'o')
%
%     % Combine hipass and lopass tapers:
%     t=sort(rand(100,1));
%     plot(t,taperfun('gausswin',t,[0 .5],4) ...
%          .* taperfun('gausswin',t,[1 .5],4),'o')
%
%    See also: TAPER, WINDOW

%     Version History:
%        Sep. 23, 2009 - initial version
%        Sep. 24, 2009 - fix no taper gives nan bug
%        Feb.  4, 2010 - added cosine taper (points to Hann), doc update
%        Feb.  5, 2012 - doc update
%        May  29, 2012 - t=limit(2) goes to taper function now
%        May  30, 2012 - flipped lim flips taper, refix nan bug, checks
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  30, 2012 at 02:20 GMT

% todo:

% check nargin
error(nargchk(3,4,nargin));

% simple checks
if(~ischar(type) || ~isvector(type))
    error('seizmo:taperfun:badInput',...
        'TYPE must be a string!');
elseif(~isnumeric(pts))
    error('seizmo:taperfun:badInput',...
        'PTS must be numeric!');
elseif(~isnumeric(limits) || ~isequal(size(limits),[1 2]))
    error('seizmo:taperfun:badInput',...
        'LIM must be a numeric vector as [lo_bound hi_bound]!');
end

% hp or lp taper
if(diff(limits)>=0)
    % hi pass
    z=pts<limits(1);
    u=pts>limits(2);
else
    % lo pass
    z=pts>limits(1);
    u=pts<limits(2);
end

% pts to be zero
if(any(z)); pts(z)=0; end

% untapered pts
if(any(u)); pts(u)=1; end

% tapered
t=~u & ~z;
if(any(t))
    func=str2func(['taper_' type]);
    if(diff(limits)==0) % no taper if width=0
        pts(t)=1;
    else
        pts(t)=(pts(t)-limits(1))/(limits(2)-limits(1));
    end
    pts(t)=func(pts(t),varargin{:});
end

end

function [pts]=taper_barthannwin(pts,varargin)
pts=0.5*pts-0.5;
pts=0.62-0.48*abs(pts)+0.38*cos(2*pi*pts);
end

function [pts]=taper_bartlett(pts,varargin)
% same as triangular
pts=taper_triang(pts);
end

function [pts]=taper_blackman(pts,varargin)
pts=0.42-0.5*cos(pi*pts)+0.08*cos(2*pi*pts);
end

function [pts]=taper_blackmanharris(pts,varargin)
pts=0.35875-0.48829*cos(pi*pts) ...
    +0.14128*cos(2*pi*pts)-0.01168*cos(3*pi*pts);
end

function [pts]=taper_bohmanwin(pts,varargin)
pts=pts.*cos(pi*(1-pts))+1/pi*sin(pi*(1-pts));
end

function [pts]=taper_chebwin(pts,r)
if(nargin==1 || isempty(r) || isnan(r)); r=100; end
% pass to signal toolbox and interpolate
w=chebwin(201,r); w=w(1:101);
pts=interp1(0:0.01:1,w,pts);
end

function [pts]=taper_cosine(pts,varargin)
pts=taper_hann(pts);
end

function [pts]=taper_flattopwin(pts,varargin)
pts=0.21557895-0.41663158*cos(pi*pts)+0.277263158*cos(2*pi*pts) ...
    -0.083578947*cos(3*pi*pts)+0.006947368*cos(4*pi*pts);
end

function [pts]=taper_gausswin(pts,a)
if(nargin==1 || isempty(a) || isnan(a)); a=2.5; end
pts=exp(-1/2*(a*(pts-1)).^2);
end

function [pts]=taper_hamming(pts,varargin)
pts=0.54-0.46*cos(pi*pts);
end

function [pts]=taper_hann(pts,varargin)
pts=(1-cos(pi*pts))/2;
end

function [pts]=taper_kaiser(pts,beta)
if(nargin==1 || isempty(beta) || isnan(beta)); beta=0.5; end
bes=abs(besseli(0,beta));
pts=besseli(0,beta*sqrt(1-(1-pts).^2))/bes;
end

function [pts]=taper_nuttallwin(pts,varargin)
pts=0.3635819-0.4891775*cos(pi*pts) ...
    +0.1365995*cos(2*pi*pts)-0.0106411*cos(3*pi*pts);
end

function [pts]=taper_parzenwin(pts,varargin)
t=pts<0.5;
if(any(t))
    pts(t)=2*pts(t).^3;
end
t=~t;
if(any(t))
    pts(t)=1-6*(1-pts(t)).^2.*pts(t);
end
end

function [pts]=taper_rectwin(pts,varargin)
pts(:)=1;
end

function [pts]=taper_triang(pts,varargin)
end

function [pts]=taper_tukeywin(pts,r)
if(nargin==1 || isempty(r) || isnan(r)); r=0.5; end
if(r<=0)
    pts(:)=1;
elseif(r>=1)
    pts=taper_hann(pts);
else
    t=pts<r;
    if(any(t))
        pts(t)=taper_hann(pts(t)/r);
    end
    t=~t;
    if(any(t))
        pts(t)=1;
    end
end
end
