function [psdgram]=read_ndbc_swden(varargin)
%READ_NDBC_SWDEN    Reads NDBC bouy spectral power density data
%
%    Usage:    psdgram=read_ndbc_swden(file)
%
%    Description:
%     PSDGRAM=READ_NDBC_SWDEN(FILE) reads the uncompressed National Data
%     Bouy Center (NDBC) power spectral density data in ascii file FILE.
%     The output PSDGRAM is a struct with the following fields:
%       PSDGRAM.name    - filename
%              .time    - timing of measurements (in DATENUM format)
%              .freq    - frequencies of measurements (in Hz)
%              .spectra - power spectral density in m^2/Hz
%              .units   - units of spectra
%
%    Notes:
%     - Ocean bouy power spectral density data are made available by the
%       National Data Bouy Center (NDBC) as compressed ascii files
%       (e.g., .txt.gz) at the following website:
%           http://dods.ndbc.noaa.gov/thredds/catalog/data/swden/
%       This function reads the uncompressed ascii files (see GUNZIP).
%
%    Examples:
%     % Read and plot the PSD for bouy 45001 from a
%     % series of storms in the middle of October 2011:
%     file=which('45001w2011101300-2011102300.txt');
%     psdgram=read_ndbc_swden(file);
%     plotpsdgram(psdgram);
%
%    See also: PLOTPSDGRAM

%     Version History:
%        Feb. 18, 2013 - doc update
%        Apr. 15, 2013 - add .name field (filename)
%        Aug. 28, 2013 - allow multiple files, fixed example
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 28, 2013 at 13:30 GMT

% todo:

% set filterspec appropriately
global SEIZMO
SEIZMO.ONEFILELIST.FILTERSPEC={'*.txt;*.TXT' 'TXT Files (*.txt,*.TXT)';
    '*.*' 'All Files (*.*)'};

% compile file lists
varargin=onefilelist(varargin{:});
nfiles=numel(varargin);

% error if no files
if(nfiles<1)
    error('seizmo:read_ndbc_swden:noFilesToRead','No files to read!');
end

% pre-allocating psdgram structure
psdgram(nfiles,1)=struct('name',[],'time',[],'freq',[],...
    'spectra',[],'units','m^2/Hz');

% read in and parse info
% - no checking done
for i=1:nfiles
    psdgram(i).name=[varargin(i).path varargin(i).name];
    d=importdata(psdgram(i).name,' ',1);
    psdgram(i).time=datenum([d.data(:,1:5) zeros(size(d.data,1),1)]);
    psdgram(i).freq=str2double(d.colheaders(6:end)).';
    psdgram(i).spectra=d.data(:,6:end);
end

end
