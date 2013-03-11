function [report]=chkpsdgram(psdgram)
%CHKPSDGRAM    Validate a psdgram struct
%
%    Usage:    msg=chkpsdgram(psdgram)
%
%    Description:
%     MSG=CHKPSDGRAM(PSDGRAM) checks that PSDGRAM is a struct as output by
%     READ_NDBC_SWDEN and can be plotted by PLOTPSDGRAM.  See those
%     functions for more details.  MSG is an error structure if a problem
%     is found and an empty matrix otherwise.
%
%    Notes:
%
%    Examples:
%     % Validate a power spectral density spectrogram
%     % data for bouy 45001 during the year 2011:
%     path=fileparts(which('read_ndbc_swden'));
%     s=read_ndbc_swden([path filesep '45001w2011.txt']);
%     error(chkpsdgram(s));
%
%    See also: PLOTPSDGRAM

%     Version History:
%        Feb. 27, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 27, 2013 at 13:30 GMT

% todo:

% check nargin
report=[];
error(nargchk(1,1,nargin));

% require scalar struct
if(~isstruct(psdgram) || ~isscalar(psdgram))
    report.identifier='seizmo:chkpsdgram:dataNotStruct';
    report.message='PSDGRAM must be a scalar struct!';
    return;
end

% check model has the required fields
reqfields={'time' 'freq' 'spectra' 'units'};
if(~all(ismember(reqfields,fieldnames(psdgram))))
    report.identifier='seizmo:chkpsdgram:badFields';
    report.message=[...
        sprintf('PSDGRAM must be a struct with fields:\n') ...
        sprintf('%s ',reqfields{:})];
    return;
end

% check fields
tsz=size(psdgram.time);
fsz=size(psdgram.freq);
if(~isreal(psdgram.time) || numel(tsz)>2 || tsz(2)~=1)
    report.identifier='seizmo:chkpsdgram:timeBad';
    report.message=['PSDGRAM .time field must be a real-valued ' ...
        'column vector!'];
    return;
elseif(~isreal(psdgram.freq) || numel(fsz)>2 || fsz(2)~=1)
    report.identifier='seizmo:chkpsdgram:freqBad';
    report.message=['PSDGRAM .freq field must be a real-valued ' ...
        'column vector!'];
    return;
elseif(~isreal(psdgram.spectra) ...
        || ~isequal([tsz(1) fsz(1)],size(psdgram.spectra)))
    report.identifier='seizmo:chkpsdgram:spectraBad';
    report.message=['PSDGRAM .spectra field must be a real-valued ' ...
        'NTxNF matrix!'];
    return;
elseif(~ischar(psdgram.units) || ndims(psdgram.units)>2 ...
        || size(psdgram.units,1)~=1)
    report.identifier='seizmo:chkpsdgram:unitsBad';
    report.message='PSDGRAM .units must be a string!';
    return;
end

end
