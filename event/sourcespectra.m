function [s,sa]=sourcespectra(Tr,Td,Mo,f)
%SOURCESPECTRA    Compute source spectra with rupture, rise & moment
%
%    Usage:    s=sourcespectra(tr,td,mo,f)
%              [s,sa]=sourcespectra(...)
%
%    Description:
%     S=SOURCESPECTRA(Td,Tr,Mo,F) returns the source amplitude spectra S at
%     the frequencies in F given the source rupture duration Tr, rise time
%     Td, scalar moment Mo.  Td, Tr, Mo, & F may be scalars or equal sized
%     arrays/vectors.
%
%     [S,SA]=SOURCESPECTRA(...) also returns the approximate source spectra
%     by using |sinc(x)|=1 for x<=1 & |sinc(x)|=1/x for x>1.
%
%    Notes:
%
%    Examples:
%     % Show the source spectra and an approximation for an
%     % earthquake with duration=100s, rise=10s, & Mo=1e25:
%     f=10.^(-4:.01:1);
%     [s,sa]=sourcespectra(100,10,1e25,f);
%     figure; loglog(f,s,f,sa);
%
%    See also: RADPAT, MAKE_SOURCE_TIMEFUNCTION

%     Version History:
%        Feb.  6, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  6, 2012 at 11:15 GMT

% todo:

% check nargin
error(nargchk(4,4,nargin));

% actual spectra
s=Mo.*abs(sin(f.*Tr/2)./(f.*Tr/2)).*abs(sin(f.*Td/2)./(f.*Td/2));

if(nargout<2); return; end

% expand scalars
[Tr,Td,Mo,f]=expandscalars(Tr,Td,Mo,f);

% make Tr the max & Td the min
[Tr,Td]=deal(max(Tr,Td),min(Tr,Td));

% approximate spectra
sa=nan(size(f));
in=f<2./Tr;
sa(in)=Mo(in);
in=~in & f<2./Td;
sa(in)=2*Mo(in)./(Tr(in).*f(in));
in=f>=2./Td;
sa(in)=4*Mo(in)./(Tr(in).*Td(in).*f(in).^2);

end
