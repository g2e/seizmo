function [snr]=quicksnr(data,nwin,swin)
%QUICKSNR    Quick estimation of SNR for SEIZMO records
%
%    Usage:    snr=quicksnr(data,noisewindow,signalwindow)
%
%    Description: QUICKSNR(DATA,NOISEWINDOW,SIGNALWINDOW) estimates the
%     signal to noise ratio for SEIZMO records by calculating the ratio of
%     the maximum-minimum amplitudes of two data windows that represent the
%     'noise' and the 'signal'.  NOISEWINDOW & SIGNALWINDOW need to be 2
%     element numeric arrays that specify the start and end of the window
%     relative to the records reference time.
%
%    Notes:
%
%    Examples:
%     To get SNR estimates of P (assuming times are stored in header):
%       Ptimes=getarrival(data,'P');
%       snr=quicksnr(data,Ptimes+[-100 -20],Ptimes+[-20 40])
%
%    See also: GETARRIVAL, CUT

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 23, 2008 - bug fix (was nsr)
%        Feb. 29, 2008 - SEISCHK support
%        Mar.  4, 2008 - doc update
%        Nov. 24, 2008 - doc update, history fix, input changed so that the
%                        windows are relative to the record reference time,
%                        better checks, formula changed to compare
%                        variation of values in the windows rather than
%                        just the maximums
%        Dec. 13, 2008 - allow different window for each record (whoops)
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:55 GMT

% todo:

% check nargin
msg=nargchk(3,3,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% turn off header checking
oldcheckheaderstate=get_checkheader_state;
set_checkheader_state(false);

% check windows
if(~isnumeric(nwin) || ~isnumeric(swin)...
        || size(nwin,2)~=2 || size(swin,2)~=2 ...
        || size(nwin,1)~=size(swin,1))
    error('seizmo:quicksnr:badInput',...
        'NOISEWINDOW & SIGNALWINDOW must be Nx2 numeric arrays!');
end

% snr=(max-min of signal)/(max-min of noise)
[nmax,nmin]=getheader(cut(data,nwin(:,1),nwin(:,2),'trim',false),'depmax','depmin');
[smax,smin]=getheader(cut(data,swin(:,1),swin(:,2),'trim',false),'depmax','depmin');
snr=(smax-smin)./(nmax-nmin);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end
