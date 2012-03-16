function [snr]=quicksnr(data,nwin,swin,method)
%QUICKSNR    Quick estimation of SNR for SEIZMO records
%
%    Usage:    snr=quicksnr(data,noisewindow,signalwindow)
%              snr=quicksnr(data,noisewindow,signalwindow,method)
%
%    Description:
%     QUICKSNR(DATA,NOISEWINDOW,SIGNALWINDOW) estimates the signal-to-noise
%     ratio for SEIZMO records by calculating the ratio of the maximum-to-
%     minimum amplitudes of two data windows that represent the 'noise' and
%     the 'signal'.  NOISEWINDOW & SIGNALWINDOW need to be 2 element
%     numeric arrays that specify the start and end of the window relative
%     to the records' reference times.
%
%     QUICKSNR(DATA,NOISEWINDOW,SIGNALWINDOW,METHOD) specifies the SNR
%     calculation method using METHOD.  METHOD must be one of the
%     following: 'PEAK2PEAK', 'RMS', 'ROBUSTRMS', 'PEAK2RMS', or
%     'PEAK2ROBUSTRMS'.  The default is 'PEAK2PEAK'.  The RMS methods
%     ('RMS' 'ROBUSTRMS') give the ratio of the RMS values.  The RMS is
%     de-meaned (or de-medianed in the 'ROBUSTRMS' case).  The last two
%     methods ('PEAK2RMS' 'PEAK2ROBUSTRMS') compare one-half the
%     peak-to-peak amplitude of the signal to the rms of the noise.
%
%    Notes:
%     - Some initial testing showed that the rms method typically gave
%       estimates that were 3/4th that of the peak2peak method while the
%       robustrms method was roughly 1/2.  This was in a rather bandlimited
%       case though so your results will differ.
%
%    Examples:
%     % To get SNR estimates of P (assuming times are stored in header):
%      P=findpicks(data,'P',1);
%      snr=quicksnr(data,P+[-100 -20],P+[-20 40])
%
%    See also: FINDPICKS, CUT

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
%        Jan. 29, 2010 - proper SEIZMO handling
%        Feb.  3, 2010 - versioninfo caching
%        Mar.  8, 2010 - versioninfo caching dropped
%        Mar. 18, 2010 - added METHOD argument (new methods: RMS, ROBUSTRMS
%                        are included), input check fix
%        Mar. 31, 2010 - added two more methods: PEAK2RMS, PEAK2ROBUSTRMS
%        Sep. 13, 2010 - nargchk fix
%        Mar. 15, 2012 - doc update, seizmocheck fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 15, 2012 at 18:00 GMT

% todo:

% check nargin
error(nargchk(3,4,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% attempt snr estimation
try
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);
    
    % check windows
    if(~isreal(nwin) || ~isreal(swin)...
            || size(nwin,2)~=2 || size(swin,2)~=2 ...
            || ~any(size(nwin,1)==[1 nrecs]) ...
            || ~any(size(swin,1)==[1 nrecs]))
        error('seizmo:quicksnr:badInput',...
            'NOISEWINDOW & SIGNALWINDOW must be 1x2 or Nx2 arrays!');
    end
    
    % check SNR method
    if(nargin==3 || isempty(method)); method='peak2peak'; end
    snrmethods={'peak2peak' 'rms' 'robustrms' 'peak2rms' 'peak2robustrms'};
    if(~ischar(method) || size(method,1)~=1 || ...
            ~any(strcmpi(method,snrmethods)))
        error('seizmo:quicksnr:badInput',...
            ['METHOD must be one of the following:\n' ...
            upper(sprintf('%s ',snrmethods{:}))]);
    end
    
    % detail message
    if(verbose); disp('Estimating SNR of Record(s)'); end
    
    % which method?
    switch lower(method)
        case 'peak2peak'
            % snr=(max-min of signal)/(max-min of noise)
            [nmax,nmin]=getheader(cut(data,nwin(:,1),nwin(:,2),...
                'trim',false),'depmax','depmin');
            [smax,smin]=getheader(cut(data,swin(:,1),swin(:,2),...
                'trim',false),'depmax','depmin');
            snr=(smax-smin)./(nmax-nmin);
        case 'rms'
            nrms=getvaluefun(cut(data,nwin(:,1),nwin(:,2),...
                'trim',false),@(x)sqrt(mean((x(:)-mean(x(:))).^2)));
            srms=getvaluefun(cut(data,swin(:,1),swin(:,2),...
                'trim',false),@(x)sqrt(mean((x(:)-mean(x(:))).^2)));
            snr=srms./nrms;
        case 'robustrms'
            nrms=getvaluefun(cut(data,nwin(:,1),nwin(:,2),...
                'trim',false),@(x)sqrt(median((x(:)-median(x(:))).^2)));
            srms=getvaluefun(cut(data,swin(:,1),swin(:,2),...
                'trim',false),@(x)sqrt(median((x(:)-median(x(:))).^2)));
            snr=srms./nrms;
        case 'peak2rms'
            % peak2peak amplitude of signal divided
            % by two compared to the rms of the noise
            nrms=getvaluefun(cut(data,nwin(:,1),nwin(:,2),...
                'trim',false),@(x)sqrt(mean((x(:)-mean(x(:))).^2)));
            [smax,smin]=getheader(cut(data,swin(:,1),swin(:,2),...
                'trim',false),'depmax','depmin');
            snr=(smax-smin)./(2*nrms);
        case 'peak2robustrms'
            % peak2peak amplitude of signal divided by
            % two compared to the robust rms of the noise
            nrms=getvaluefun(cut(data,nwin(:,1),nwin(:,2),...
                'trim',false),@(x)sqrt(median((x(:)-median(x(:))).^2)));
            [smax,smin]=getheader(cut(data,swin(:,1),swin(:,2),...
                'trim',false),'depmax','depmin');
            snr=(smax-smin)./(2*nrms);
    end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);

    % rethrow error
    error(lasterror);
end

end
