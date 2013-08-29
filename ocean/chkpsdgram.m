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
%     % data for bouy 45001 from October 2011:
%     file=which('45001w2011101300-2011102300.txt');
%     psdgram=read_ndbc_swden(file);
%     error(chkpsdgram(psdgram));
%
%    See also: PLOTPSDGRAM

%     Version History:
%        Feb. 27, 2013 - initial version
%        Apr. 15, 2013 - require .name field
%        Aug.  9, 2013 - allow psdgram to be an array
%        Aug. 28, 2013 - fixed example
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 28, 2013 at 13:30 GMT

% todo:

% check nargin
report=[];
error(nargchk(1,1,nargin));

% require struct
if(~isstruct(psdgram))
    report.identifier='seizmo:chkpsdgram:dataNotStruct';
    report.message='PSDGRAM must be a struct!';
    return;
end

% check model has the required fields
reqfields={'name' 'time' 'freq' 'spectra' 'units'};
if(~all(ismember(reqfields,fieldnames(psdgram))))
    report.identifier='seizmo:chkpsdgram:badFields';
    report.message=[...
        sprintf('PSDGRAM must be a struct with fields:\n') ...
        sprintf('%s ',reqfields{:})];
    return;
end

% check fields
for i=1:numel(psdgram)
    tsz=size(psdgram(i).time);
    fsz=size(psdgram(i).freq);
    if(~ischar(psdgram(i).name) || ndims(psdgram(i).name)~=2 ...
            || size(psdgram(i).name,1)~=1)
        report.identifier='seizmo:chkpsdgram:nameBad';
        report.message='PSDGRAM .name field must be a string!';
    elseif(~isreal(psdgram(i).time) || numel(tsz)>2 || tsz(2)~=1)
        report.identifier='seizmo:chkpsdgram:timeBad';
        report.message=['PSDGRAM .time field must be a real-valued ' ...
            'column vector!'];
        return;
    elseif(~isreal(psdgram(i).freq) || numel(fsz)>2 || fsz(2)~=1)
        report.identifier='seizmo:chkpsdgram:freqBad';
        report.message=['PSDGRAM .freq field must be a real-valued ' ...
            'column vector!'];
        return;
    elseif(~isreal(psdgram(i).spectra) ...
            || ~isequal([tsz(1) fsz(1)],size(psdgram(i).spectra)))
        report.identifier='seizmo:chkpsdgram:spectraBad';
        report.message=['PSDGRAM .spectra field must be a real-valued ' ...
            'NTxNF matrix!'];
        return;
    elseif(~ischar(psdgram(i).units) || ndims(psdgram(i).units)>2 ...
            || size(psdgram(i).units,1)~=1)
        report.identifier='seizmo:chkpsdgram:unitsBad';
        report.message='PSDGRAM .units must be a string!';
        return;
    end
end

end
