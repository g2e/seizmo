function [b]=sacsize(v,t,n,e,c)
%BYTESIZE    Returns size of SAC file in bytes
%
%    Usage:  bytes=sacsize(nvhdr,iftype,npts,leven)
%            bytes=sacsize(nvhdr,iftype,npts,leven,ncmp)

% check arguments
error(nargchk(4,5,nargin));

% header info
h=sachi(v);

% account for header size
b=h.data.startbyte;

% account for uneven data
m=0; if(e==h.false); m=1; end

% figure rest by type
if(t==h.enum(1).val.itime)
    m=m+1;
elseif(t==h.enum(1).val.irlim)
    m=m+2;
elseif(t==h.enum(1).val.iamph)
    m=m+2;
elseif(t==h.enum(1).val.ixy)
    m=m+1;
elseif(t==h.enum(1).val.ixyz)
    m=m+1;
elseif(nargin==5)
    if(c<1)
        warning('SAClab:ncmpBad','ncmp<1 not allowed for multi-component timeseries!');
        b=0; m=0;
    else
        m=m+c;
    end
end

% total
b=b+m*n*h.data.bytesize;

end