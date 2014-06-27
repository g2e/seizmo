function [f,c,a]=noise_bessel_fit(data,ctest,atest)
%NOISE_BESSEL_FIT    Fit noise data with Bessel functions
%
%    Usage:    [f,c]=noise_bessel_fit(data,ctest)
%              [f,c,a]=noise_bessel_fit(data,ctest,atest)
%
%    Description:
%     [F,C]=NOISE_BESSEL_FIT(DATA,CTEST)
%
%     [F,C,A]=NOISE_BESSEL_FIT(DATA,CTEST,ATEST)
%
%    Notes:
%
%    Examples:
%     % 
%
%    See also: NOISE_STACK_DELAZ, NOISE_STACK, STACK2STACK, NOISE_PROCESS,
%              NOISE_SETUP

%     Version History:
%        June  9, 2014 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  9, 2014 at 11:15 GMT

% todo:
% - coherency isn't
%   - cross correlation must be XY'/sqrt(XX'YY')
%     - add a whiten/coherency option?
%     - how about noise_coherency
%   - stacking must be in the frequency domain...
%     - do we need a stacking style option?
% - defaults for ctest & atest?
%   - 2-6km/s stepping at .005km/s
%   - 10^(-5:.01:-2) /km
%   - 2 outputs = no attenuation search
% - convert data to real & imag
%   - 2*delta*fft ?
% - loop over frequencies
%   - loop over ctest
%     - calculate besselj at observed distances
%     - get L1 misfit
%   - global misfit minimum gives c @ current frequency
%   - loop over atest
%     - calculate besselj at observed distances given c
%     - get logarithm of envelope for both theory and data
%       - envelope of a sparse signal
%         - interpolate missing points
%     - get L1 misfit
%       - only for theoretical values >10^-2
%   - global misfit minimum gives a @ current frequency

end
