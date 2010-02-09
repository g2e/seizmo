function [phvel]=prem_dispersion(freq)
%PREM_DISPERSION    Returns fundamental mode PREM dispersion
%
%    Usage:    phvel=prem_dispersion(freq)
%
%    Description: PHVEL=PREM_DISPERSION(FREQ) returns the Rayleigh wave
%     phase velocities PHVEL at the frequencies specified in FREQ.
%
%    Notes:
%     - At present this function is quite crappy.
%
%    Examples:
%     Plot up a PREM dispersion curve:
%     freq=0.005:0.001:0.1;
%     phvel=prem_dispersion(freq);
%     plot(freq,phvel)
%
%    See also: PREM

%     Version History:
%        Feb.  5, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  5, 2010 at 20:00 GMT

% todo
% - get better dispersion values using Saito code
% - Love waves
% - group velocities
% - quality factor

% period vs prem phase velocity
prem=[ % from a random paper
 15.0000    3.8600
 20.0000    3.9650
 25.0000    3.9900
 30.0000    4.0120
 35.0000    4.0220
 40.0000    4.0300
 45.0000    4.0350
 50.0000    4.0400
 55.0000    4.0450
 60.0000    4.0500
 65.0000    4.0580
 70.0000    4.0670
 75.0000    4.0770
 80.0000    4.0880
 85.0000    4.1010
 90.0000    4.1150
 95.0000    4.1330
100.0000    4.1500
105.0000    4.1660
110.0000    4.1820
115.0000    4.1980
120.0000    4.2140
125.0000    4.2300
130.0000    4.2460
135.0000    4.2620
140.0000    4.2780
145.0000    4.2940
150.0000    4.3100
155.0000    4.3260
160.0000    4.3420
165.0000    4.3580
170.0000    4.3740
175.0000    4.3900
180.0000    4.4060
185.0000    4.4220];

% linear interpolation
phvel=interp1(prem(:,1),prem(:,2),1./freq,[],'extrap');

end
