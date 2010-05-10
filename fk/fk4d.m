function [varargout]=fk4d(data,width,overlap,smax,spts,frng,polar,center)
%FK4D    Returns a map of energy in frequency-wavenumber-time space
%
%    Usage:    s4d=fk4d(data,width,overlap,smax,spts,frng)
%              s4d=fk4d(data,width,overlap,smax,spts,frng,polar)
%              s4d=fk4d(data,width,overlap,smax,spts,frng,polar,center)
%
%    Description: S4D=FK4D(DATA,WIDTH,OVERLAP,SMAX,SPTS,FRNG) calculates
%     the energy passing through an array in frequency-wavenumber-time
%     space.  Actually, to allow for easier interpretation between
%     frequencies, the energy is mapped into frequency-slowness-time space.
%     The addition of the time dimension (compared to FKVOLUME) is made
%     possible by taking windows that are a fraction of the input records
%     in SEIZMO struct DATA.  The array info is derived from the metadata
%     stored in the SEIZMO struct DATA, so make sure station location and
%     timing fields are set!  WIDTH is the width of the time windows in
%     percent of the total time extent of the records & OVERLAP is the
%     percent of the overlap between each window.  The default WIDTH is
%     1% and the default OVERLAP is 0%.  The range of the slowness space is
%     specified by SMAX (in s/deg) and extends from -SMAX to SMAX for both
%     East/West and North/South directions.  SPTS controls the number of
%     slowness points for both directions (SPTSxSPTS grid).  FRNG gives the
%     frequency range as [FREQLOW FREQHIGH] in Hz.  S4D is a struct array
%     whose elements are essentially time frames of the frequency-slowness
%     space and contain relevant info as well as the frequency-slowness
%     data.  The struct layout is:
%          .response - frequency-slowness array response
%          .nsta     - number of stations utilized in making map
%          .stla     - station latitudes
%          .stlo     - station longitudes
%          .stel     - station elevations (surface)
%          .stdp     - station depths (from surface)
%          .butc     - UTC start time of data
%          .eutc     - UTC end time of data
%          .npts     - number of time points
%          .delta    - number of seconds between each time point
%          .x        - east/west slowness or azimuth values
%          .y        - north/south or radial slowness values
%          .z        - frequency values
%          .polar    - true if slowness is sampled in polar coordinates 
%          .center   - array center or method
%          .normdb   - what 0dB actually corresponds to
%          .volume   - true if frequency-slowness volume (false for FKMAP)
%
%     S4D=FK4D(DATA,WIDTH,OVERLAP,SMAX,SPTS,FRNG,POLAR) specifies if the
%     slowness space is sampled regularyly in cartesian or polar
%     coordinates.  Polar coords are useful for slicing the volume by
%     azimuth (pie slice) or slowness magnitude (rings).  Cartesian coords
%     (the default) samples the slowness space regularly in the East/West &
%     North/South directions and so exhibits less distortion of the
%     slowness space.
%
%     S4D=FK4D(DATA,WIDTH,OVERLAP,SMAX,SPTS,FRNG,POLAR,CENTER) defines the
%     array center.  CENTER may be [LAT LON], 'center', 'coarray', or
%     'full'.  The default is 'coarray'.  The 'center' option finds the
%     center position of the array by averaging the station positions
%     (using ARRAYCENTER).  Both 'coarray' and 'full' are essentially
%     centerless methods using the relative positioning between every
%     possible pairing of stations in the array.  The 'full' method
%     includes redundant and same station pairings (and will always give
%     poorer results compared to 'coarray').
%
%    Notes:
%     - Records in DATA must have equal number of points, equal sample
%       spacing, the same start time (in absolute time), and be evenly
%       spaced time series records.  Use functions SYNCHRONIZE, SYNCRATES,
%       & INTERPOLATE to get the timing/sampling the same.
%
%    Examples:
%     Get frequency-slowness-time data for an array at 20-50s periods:
%      s4d=fk4d(data,1,75,50,201,[1/50 1/20]);
%
%    See also: PLAYFKVOLUME, FKMAP, FKARF,
%              SNYQUIST, KXY2SLOWBAZ, SLOWBAZ2KXY

%     Version History:
%        May   9, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   9, 2010 at 15:30 GMT

% todo:
% - require 6 to 8 inputs
% - check the dataset (needs fields defined)
% - check inputs
% - get earliest and latest times
% - interpolate window based on those times and width and overlap
% - do we set data outside window to zero or do we throw it out?
%   - i think dropping would be more accurate (if not the same)
% - then just pass to fkvolume
% - the end

end
