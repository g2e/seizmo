function [psdgram]=read_ndbc_swden(file)
%READ_NDBC_SWDEN    Reads NDBC bouy spectral power density data
%
%    Usage:    psdgram=read_ndbc_swden(file)
%
%    Description:
%     PSDGRAM=READ_NDBC_SWDEN(FILE) reads the uncompressed National Data
%     Bouy Center (NDBC) power spectral density data in ascii file FILE.
%     The output PSDGRAM is a struct with the following fields:
%       PSDGRAM.time    - timing of measurements (in DATENUM format)
%              .freq    - frequencies of measurements (in Hz)
%              .spectra - power spectral density in m^2/Hz
%              .units   - units of spectra
%
%    Notes:
%     - Ocean bouy power spectral density data are made available by the
%       National Data Bouy Center (NDBC) as gunzipped ascii files.  This
%       reads
%
%    Examples:
%     % Read and plot PSD for bouy 45001 in 2011:
%     psdgram=read_ndbc_swden('45001w2011.txt');
%     plotpsdgram(psdgram);
%
%    See also: PLOTPSDGRAM

%     Version History:
%        Feb. 18, 2013 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 18, 2013 at 13:30 GMT

% todo:

error(nargchk(0,1,nargin));
if(~nargin); file=[]; end
d=importdata(file,' ',1);
psdgram.time=datenum([d.data(:,1:5) zeros(size(d.data,1),1)]);
psdgram.freq=str2double(d.colheaders(6:end)).';
psdgram.spectra=d.data(:,6:end);
psdgram.units='m^2/Hz';

end
