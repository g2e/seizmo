function [s,sa]=sourcespectra(Tr,Td,Mo,f)
%SOURCESPECTRA    Compute source spectra with rupture, rise & moment
%
%    Usage:    s=sourcespectra(tr,td,mo,f)
%              [s,sa]=sourcespectra(...)
%
%    Description:
%     S=SOURCESPECTRA(Tr,Td,Mo,F) returns the source amplitude spectra S of
%     an earthquake at the frequencies in F given the source rupture
%     duration Tr, rise time Td, and scalar moment Mo.  The spectra is
%     derived from a simple model of the source time function as the
%     convolution of 2 unit area boxcar functions of length Tr & Td and
%     scaling the result by Mo.  This is easily implemented in the
%     frequency domain as the sinc function is the Fourier transform pair
%     of a Boxcar function.  Multiplication of the moment Mo by the sinc
%     function for rupture time Tr as well as the sinc function for rise
%     time Td gives the amplitude spectra.  Td, Tr, Mo, & F may be scalars
%     or equal sized arrays/vectors.
%
%     [S,SA]=SOURCESPECTRA(...) also returns an approximate source spectra
%     by using |sinc(f)|=1 for f<=2/Tx & |sinc(f)|=1/f for f>2/Tx instead
%     of the sinc functions of Td & Tr.  This allows clearly seeing the
%     corner frequencies at 2/Tr & 2/Td.
%
%    Notes:
%     - For an explanation of Tr & Td see Stein & Wysession 2002 sections
%       4.3.2 & 4.6.2.
%     - Formulas from Stein & Wysession 2002 page 267 equations 7 & 9.
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
%        Feb. 18, 2014 - doc update, input checking
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 18, 2014 at 11:15 GMT

% todo:

% check nargin
error(nargchk(4,4,nargin));

% check inputs
if(~isequalsizeorscalar(Tr,Td,Mo,f))
    error('seizmo:sourcespectra:badInput',...
        'Inputs must be equal sized or scalar!');
end

% actual spectra
s=Mo.*abs(sin(f.*Tr/2)./(f.*Tr/2)).*abs(sin(f.*Td/2)./(f.*Td/2));

% do we want an approximated source spectra?
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
