function [order,wn] = buttord2(wp,ws,rp,rs,opt)
%BUTTORD2 Butterworth filter order selection. (Honors passband corners)
%   [N, Wn] = BUTTORD2(Wp, Ws, Rp, Rs) returns the order N of the lowest 
%   order digital Butterworth filter that loses no more than Rp dB in 
%   the passband and has at least Rs dB of attenuation in the stopband.  
%   Wp and Ws are the passband and stopband edge frequencies, normalized 
%   from 0 to 1 (where 1 corresponds to pi radians/sample). For example,
%       Lowpass:    Wp = .1,      Ws = .2
%       Highpass:   Wp = .2,      Ws = .1
%       Bandpass:   Wp = [.2 .7], Ws = [.1 .8]
%       Bandstop:   Wp = [.1 .8], Ws = [.2 .7]
%   BUTTORD2 also returns Wn, the Butterworth natural frequency (or, 
%   the "3 dB frequency") to use with BUTTER to achieve the specifications. 
%
%   [N, Wn] = BUTTORD2(Wp, Ws, Rp, Rs, 's') does the computation for an 
%   analog filter, in which case Wp and Ws are in radians/second.
%
%   When Rp is chosen as 3 dB, the Wn in BUTTER is equal to Wp in BUTTORD2.
%   This is not the case for BUTTORD which honors stopband corner
%   positioning over passband corner positioning.
%
%   See also BUTTORD, BUTTER, CHEB1ORD, CHEB2ORD, ELLIPORD.

%   Author(s): L. Shure, 6-9-88
%              T. Krauss, 11-13-92, revised
%              G. Euler, 2010-02-02, require exactly rp at wp
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.10.4.2 $  $Date: 2004/12/26 22:15:26 $

%   Reference(s):
%     [1] Rabiner and Gold, p 241.

error(nargchk(4,5,nargin));
error(nargoutchk(0,2,nargout));

if nargin == 4
	opt = 'z';
elseif nargin == 5
	if ~strcmp(opt,'z') && ~strcmp(opt,'s')
		error('Invalid option for final argument.');
	end
end

msg=freqchk(wp,ws,opt);
error(msg);

ftype = 2*(length(wp) - 1);
if wp(1) < ws(1)
	ftype = ftype + 1;	% low (1) or reject (3)
else
	ftype = ftype + 2;	% high (2) or pass (4)
end

% first, prewarp frequencies from digital (unit circle) to analog (imag. axis):
if strcmp(opt,'z')	% digital
	WP=tan(pi*wp/2);
	WS=tan(pi*ws/2);
else  % don't have to if analog already
	WP=wp;
	WS=ws;
end
%note - on old systems that are NOT case sensitive, this will still work OK

% next, transform to low pass prototype with passband edge of 1 and stopband
% edges determined by the following: (see Rabiner and Gold, p.258)
if ftype == 1	% low
	WA=WS/WP;
elseif ftype == 2	% high
	WA=WP/WS;
elseif ftype == 3	% stop
	WA=(WS*(WP(1)-WP(2)))./(WS.^2 - WP(1)*WP(2));
elseif ftype == 4	% pass
	WA=(WS.^2 - WP(1)*WP(2))./(WS*(WP(1)-WP(2)));
end


% find the minimum order b'worth filter to meet the more demanding spec:
WA=min(abs(WA));
order = ceil( log10( (10 .^ (0.1*abs(rs)) - 1)./ ...
	(10 .^ (0.1*abs(rp)) - 1) ) / (2*log10(WA)) );

% next find the butterworth natural frequency W0 (or, the "3dB frequency")
% to give exactly rp dB at 1.  W0 will be between 1 and WA:
W0 = 1 / ( (10^(.1*abs(rp)) - 1)^(1/(2*(abs(order)))));

% now convert this frequency back from lowpass prototype 
% to the original analog filter:
if ftype == 1	% low
	WN=W0*WP;
elseif ftype == 2	% high
	WN=WP/W0;
elseif ftype == 3	% stop
	WN(1) = ( (WP(2)-WP(1)) + sqrt((WP(2)-WP(1))^2 + ...
		4*W0.^2*WP(1)*WP(2)))./(2*W0);
	WN(2) = ( (WP(2)-WP(1)) - sqrt((WP(2)-WP(1))^2 + ...
		4*W0.^2*WP(1)*WP(2)))./(2*W0);
	WN=sort(abs(WN)); 
elseif ftype == 4	% pass
	W0=[-W0 W0];  % need both left and right 3dB frequencies
	WN= -W0*(WP(2)-WP(1))/2 + sqrt( W0.^2/4*(WP(2)-WP(1))^2 + WP(1)*WP(2) );
	WN=sort(abs(WN)); 
end

% finally, transform frequencies from analog to digital if necessary:
if strcmp(opt,'z')	% digital
	wn=(2/pi)*atan(WN);  % bilinear transform
else
	wn=WN;
end

function errmsg=freqchk(wp,ws,opt)
%FREQCHK Parameter checking for BUTTORD, CHEB1ORD, CHEB2ORD, ELLIPORD.
%
%   MSG=FREQCHK(WP,WS,OPT) checks for the correct passband (WP) and 
%   stopband (WS) frequency specifications and returns a diagnostic 
%   message (MSG). If the parameters are specified correctly MSG is
%   empty. OPT is a string indicating digital or analog filter design.
%
%   The frequency parameters are checked to be of the same length.
%   For bandpass and bandstop filters, the frequency bands are checked
%   for increasing values and no overlaps.
% 
%   For digital filters, WP and WS are checked to be in (0,1). 
%   For analog filters, WP and WS are checked to be non-negative values.
   
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 2002/04/15 01:08:57 $

errmsg='';

% Check for correct lengths
if length(wp)~=length(ws),
   errmsg = 'The frequency vectors must both be the same length.';
   return
end

% Check for allowed interval values 
if strcmp(opt,'z'),
   if any(wp<=0) || any(wp>=1) || any(ws<=0) || any(ws>=1),
      errmsg = 'The cutoff frequencies must be within the interval of (0,1).';
      return    
   end
else % Analog filter design
   if any(wp<=0) || any(ws<=0),
      errmsg = 'The cutoff frequencies must be non-negative for analog filters.';
      return    
   end
end

% For Band specifications
if length(wp)==2,
   % Check for frequencies to be in increasing order
   if (wp(1)>=wp(2)) || (ws(1)>=ws(2)),
      errmsg = 'The cutoff frequencies should be in increasing order.';
      return
   end    
   % Check for passband and stopband frequency overlaps  
   if ~( ((wp(1)<ws(1)) && (wp(2)>ws(2))) || ((wp(1)>ws(1)) && (wp(2)<ws(2))) ),
      errmsg = 'The passband and stopband cutoff frequencies should not overlap.'; 
      return
   end
end

% [EOF] freqchk.m
