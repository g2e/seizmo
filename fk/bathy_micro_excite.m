function [varargout]=bathy_micro_excite(z,f,vs)
%BATHY_MICRO_EXCITE    Bathymetry 2ndary microseism excitation coefficients
%
%    Usage:    c=bathy_micro_excite(z,f,vs)
%              [c1,c2,c3,c4]=bathy_micro_excite(z,f,vs)
%
%    Description:
%     C=BATHY_MICRO_EXCITE(Z,F,Vs) returns the summed squares of the
%     bathymetry secondary microseism excitation coefficients using the
%     values listed in Table 1 of Longuet-Higgins 1950 (pg. 30).  This
%     modulates the wave-wave interaction intensity of oceanic gravity
%     waves to give the microseism energy (Kedar et al 2008).  Z is the
%     bathymetry and may be an array (C is the same size as Z).  F is the
%     frequency (not the angular frequency) and Vs is the shear velocity of
%     the sea bed in meters per second.  The default value for F is 1/(2pi)
%     and Vs is 2800 m/s.
%
%     [C1,C2,C3,C4]=BATHY_MICRO_EXCITE(Z,F,VS) returns the individual
%     components from C1 to C4.  These can be used to calculate the former
%     usage output C=C1.^2+C2.^2+C3.^2+C4.^2.
%
%    Notes:
%     - References:
%        Longuet-Higgins 1950, Theory of the Origin of Microseisms
%        Kedar et al 2008, The origin of Deep Ocean Microseisms in the
%         North Atlantic Ocean
%
%    Examples:
%     % Crust2.0 excitation:
%     [lon,lat]=meshgrid(-179:2:179,89:-2:-89);
%     c2elev=getc2elev(lat,lon);
%     c2elev(c2elev>0)=0; % mask out land
%     c=bathy_micro_excite(-c2elev);
%     figure;
%     imagesc(lon(1,:),lat(:,1),c);
%     set(gca,'ydir','normal');
%
%    See also: SLOWNESS2DEG

%     Version History:
%        Sep. 22, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 22, 2010 at 10:35 GMT

% todo

% check nargin
error(nargchk(1,3,nargin));

% defaults
if(nargin<2 || isempty(f)); f=1/(2*pi); end
if(nargin<3 || isempty(vs)); vs=2800; end

% check inputs
if(~isreal(z) || any(z(:)<0))
    error('seizmo:bathy_micro_excite:badInput',...
        'Z must be the bathymetry (positive reals in meters)!');
elseif(~isreal(f) || any(f(:)<0))
    error('seizmo:bathy_micro_excite:badInput',...
        'F must be the frequency (positive reals in Hz)!');
elseif(~isreal(vs) || any(vs(:)<100))
    error('seizmo:bathy_micro_excite:badInput',...
        'Vs must be the sea bed shear velocity (positive reals in m/s)!');
end

% bathymetry modulated excitation coefficients
% Column 1 = ang freq * bathymetry / shear velocity of sea floor
% Column 2 = coefficient of excitation
c1=[0.00 .191
    0.10 .206
    0.48 .368
    0.63 .565
    0.72 .728
    0.77 .837
    0.82 .894
    0.85 .908
    0.89 .890
    0.92 .857
    0.99 .759
    1.06 .649
    1.13 .542
    1.22 .444
    1.33 .355
    1.48 .276
    1.68 .205
    1.99 .139
    2.59 .078
    3.23 .049
    3.87 .034
    4.87 .021
    5.31 .017
    5.92 .014];
c2=[1.01 .000
    1.03 .038
    1.04 .076
    1.06 .108
    1.09 .141
    1.14 .168
    1.20 .178
    1.28 .180
    1.36 .177
    1.46 .173
    1.58 .170
    1.72 .172
    1.86 .180
    1.98 .194
    2.39 .318
    2.58 .418
    2.70 .454
    2.80 .448
    2.89 .421
    2.97 .386
    3.06 .351
    3.14 .316
    3.33 .256
    3.54 .206
    3.79 .165
    4.09 .131
    4.47 .101
    4.99 .076
    5.73 .054
    6.96 .034];
c3=[2.83 .000
    2.84 .036
    2.86 .070
    2.88 .098
    2.91 .126
    2.97 .151
    3.04 .163
    3.11 .166
    3.21 .165
    3.32 .163
    3.44 .162
    3.58 .164
    3.73 .171
    3.86 .184
    4.30 .280
    4.53 .331
    4.69 .330
    4.83 .305
    4.96 .275
    5.09 .245
    5.22 .218
    5.36 .194
    5.67 .154
    6.02 .123];
c4=[4.64 .000
    4.65 .031
    4.67 .065
    4.69 .090
    4.73 .115
    4.80 .138
    4.87 .150
    4.95 .155
    5.05 .156
    5.17 .155
    5.30 .154
    5.44 .157
    5.60 .164
    5.73 .175
    6.21 .250];

% normalized
x=2*pi*f.*z./vs;

% spline interpolation
sc1=interp1(c1(:,1),c1(:,2),x,'spline',0);
sc2=interp1(c2(:,1),c2(:,2),x,'spline',0);
sc3=interp1(c3(:,1),c3(:,2),x,'spline',0);
sc4=interp1(c4(:,1),c4(:,2),x,'spline',0);

% nargout
if(nargout>1)
    varargout={sc1 sc2 sc3 sc4};
else
    varargout={sc1.^2+sc2.^2+sc3.^2+sc4.^2};
end

end
