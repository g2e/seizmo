function [cg,lg,pg]=mcxc(x,varargin)
%MCXC    Multi-channel cross correlation with built-in fast peak picker
%
%    Usage:   [cg,lg]=mcxc(x,y)
%             [cg,lg]=mcxc(x)
%             [cg,lg,pg]=mcxc(...,'npeaks',npeaks,...)
%             [cg,lg,pg]=mcxc(...,'spacing',spacing,...)
%             [cg,lg,pg]=mcxc(...,'lags',lags,...)
%             [cg,lg,pg]=mcxc(...,'adjacent',adjacent,...)
%             [cg,lg,pg]=mcxc(...,'normxc',normxc,...)
%             [cg,lg,pg]=mcxc(...,'absxc',absxc,...)
%             [cg,lg,pg]=mcxc(...,'vectorxc',vectorxc,...)
%             [cg,lg,pg]=mcxc(...,'pow2pad',pow2pad,...)
%             [cg,lg,pg]=mcxc(...,'verbose',verbosity,...)
%
%    Description:
%     [CG,LG]=MCXC(X,Y) performs multi-channel cross correlation of records
%     distributed in the columns of X against records distributed in the
%     columns of Y returning the correlograms in CG such that the row and
%     column correpond to the columns of the two records being cross
%     correlated in X and Y.  The entire correlogram is stored down the
%     third dimension of CG such that CG(3,7,:) would give the correlogram
%     of cross correlating record 3 in X against record 7 in Y.  LG is a
%     vector giving the corresponding lags.
%
%     [CG,LG]=MCXC(X) performs multi-channel cross-correlation of every
%     unique pairing for the records distributed in the columns of X.
%     Autocorrelations are not included.  CG contains the correlograms
%     distributed such that each has its own row and they are stored down
%     the third dimension so CG(1,1,:) gives the correlogram for records 1
%     and 2, CG(2,1,:) gives the correlogram for records 1 and 3, and so
%     on. LG is a vector giving the corresponding lags.
%
%     [CG,LG,PG]=MCXC(...,'npeaks',NPEAKS,...) controls the number of peaks
%     to pick in each correlogram.  By default NPEAKS is set to 0 which
%     causes MCXC to return the entire correlogram.  Setting NPEAKS to 1 or
%     more will turn on the peak picker (which picks the maximum valued
%     point in the correlogram within the constraints set through the
%     options and so may not always pick a 'peak') and only returns the
%     correlogram value, lag, and polarity of NPEAKS peaks in CG, LG, and
%     PG respectively.  Multiple peaks will be stored down the third
%     dimension of the returned matrices such that CG(:,:,N), LG(:,:,N),
%     PG(:,:,N) give info about the Nth highest peak of every correlogram.
%
%     [CG,LG,PG]=MCXC(...,'spacing',SPACING,...) controls the minimum
%     spacing between returned peaks.  By default SPACING is set to 1 which
%     requires that peaks be at least 1 sample interval apart.  Setting
%     SPACING==0 will cause the peak picker to always return the same peak
%     for all the peaks.  Setting SPACING too high can cause the peak
%     picker to return a peak with value=nan, lag=?, polarity=nan when no
%     points are left satisfying the SPACING requirement.  SPACING will be 
%     ignored unless NPEAKS is set >0.
%
%     [CG,LG,PG]=MCXC(...,'lags',LAGS,...) controls the lag range that the
%     peak picker is allowed to search.  By default LAGS is set to the
%     maximum possible range.  LAGS can be a one (symmetric window) or two
%     (asymmetric window) element vector indicating the range of lags (in 
%     samples!) to search.
%
%     [CG,LG,PG]=MCXC(...,'adjacent',ADJACENT,...) controls the number of
%     adjacent points to be returned in addition to each pick picked.  By
%     default ADJACENT is set to 0 which requires that all points within 0
%     sample intervals of each peak be returned.  Setting ADJACENT to 2
%     will return the peak and all points at and within 2 sample intervals
%     bringing the total to 5 points for each peak.  The peaks and their
%     adjacent points info are stored down the 4th dimension of CG, LG, and
%     PG and are kept in their original order such that the peak is always
%     in the middle.  For example, with ADJACENT set to 3, CG(:,:,1,4)
%     would give the highest peak cross correlation value for each 
%     correlogram (note that floor(ADJACENT)+1 will always give the 4th
%     dimension index of the peak values).  Adjacent points may be set to 
%     value=nan, lag=?, polarity=nan when the point comes within SPACING of
%     another peak or if the point falls outside the allowed lag range.  
%     ADJACENT will be ignored unless NPEAKS is set >0.
%
%     [CG,LG,PG]=MCXC(...,'normxc',NORMXC,...) controls whether the cross
%     correlations are to be normalized (using the autocorrelations) or
%     not.  By default NORMXC is set to TRUE.  Correlation values for a
%     normalized correlogram are in the range from -1 to 1, with 1 being a
%     perfect correlation and -1 being a perfect anticorrelation.
%
%     [CG,LG,PG]=MCXC(...,'absxc',ABSXC,...) controls whether the peak
%     picker looks at absolute peaks or just positive peaks.  By default
%     ABSXC is set to TRUE which causes the peak picker to first take the
%     absolute value of the correlograms before picking the highest peaks.
%     Polarities are returned in the PG matrix.  Setting ABSXC to FALSE 
%     will cause the peak picker to search for the highest peak in the 
%     unaltered correlograms and so PG will always just contain 1's.  ABSXC
%     is ignored if NPEAKS==0.
%
%     [CG,LG,PG]=MCXC(...,'vectorxc',VECTORXC,...) controls which
%     subfunctions are used in the computation: those that use for loops or
%     those that have been 'vectorized'.  By default VECTORXC is set to
%     false which means that computations are done using the for looped
%     versions.  This is mainly for development or benchmarking.  Currently
%     the vectorized versions are SLOWER due to Matlab's JIT compiler
%     preferring small looped operations over chunky vectorized operations.
%     Vectorized versions also have a significantly larger memory footprint
%     which may cause issues with some systems.  Changing VECTORXC is 
%     discouraged for these reasons.
%
%     [CG,LG,PG]=MCXC(...,'pow2pad',POW2PAD,...) controls the length of
%     zero padding for spectral operations and so controls the frequency 
%     resolution.  FFT length is set according to:
%         fftlength=2^(nextpow2(max([size(X,1) size(Y,1)]))+POW2PAD)
%     POW2PAD may be any positive integer where higher integers give higher
%     spectral resolution for the correlation.  By default POW2PAD is set
%     to 1 and changing it is discouraged.
%
%     [CG,LG,PG]=MCXC(...,'verbose',VERBOSITY,...) switches the verbose
%     message output for MCXC.  VERBOSITY should be a logical where TRUE
%     gives verbose output and FALSE does not.  The default is TRUE.
%
%    Notes:
%     Options Summary:
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Defaults:                                Valid Range:
%      - NORMXC   == TRUE     Any value that can be evaluated as TRUE/FALSE
%      - ABSXC    == TRUE     Any value that can be evaluated as TRUE/FALSE
%      - VECTORXC == FALSE    Any value that can be evaluated as TRUE/FALSE
%      - NPEAKS   == 0                          0+ (INTEGERS)
%      - SPACING  == 1                          0+
%      - ADJACENT == 0                          0+
%      - LAGS     == [-size(X,1)+1 size(Y,1)-1] -size(X,1)+1 to size(Y,1)-1
%      - POW2PAD  == 1                          1+ (INTEGERS)
%      - VERBOSE  == TRUE     Any value that can be evaluated as TRUE/FALSE
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Output Summary (with NPEAKS >0):
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      CG,LG,PG SIZE
%        [size(X,2),size(Y,2),NPEAKS,ADJACENT] if MCXC(X,Y)
%                    or
%        [size(X,2)*(size(X,2)-1)/2,1,NPEAKS,ADJACENT] if MCXC(X)
%
%      CG: Contains the correlation values at peaks ordered from highest to
%       lowest down the third dimension and with adjacent points down the
%       4th dimension. ADJACENT+1 gives the 4th dimension index for the 
%       peaks.
%      LG: Lags corresponding to values given in CG.
%      PG: Signs corresponding to values in CG.  Useful if ABSXC is TRUE.
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Examples:
%     % The single dataset option is good for cases where correlations
%     % between all records are needed.  For instance, in a seismic noise
%     % study.  This will produce correlograms for all pairs of records
%     % in a dataset (Remember individual records should extend down the
%     % columns of X!):
%     [correlograms,lags]=mcxc(X);
%
%    See also: XCORR, FFTFILT (Signal Processing Toolbox)

%    Seismology Specific Notes:
%     The idea behind MCXC is that it is a more robust way of determining
%     the relative delays of phases between stations than using a 'master'
%     station system because it reduces the adverse effects of noise and
%     waveform distortion.  It also doesn't require the user to select a
%     master station in the first place (which is a crucial part of that
%     method) - a boon for automation.  Probably the largest benefit to the
%     mcxc method is that the results can be incorporated into a cluster
%     analysis to quantitatively subset the data as other later analytical 
%     methods may require, or to automatically remove noisy records, or to
%     search for oddball recordings (this is a great way to search for 
%     stations with improper instrument response).  The main issue is that
%     the mcxc method requires significantly more computation time and the
%     results need to be pushed through an inversion to get the delay times
%     of phase arrivals between stations.
%
%     This particular version implements a multi-peak picker, which is
%     useful for dealing with noisy datasets.  Use the NPEAKS, ADJACENT and 
%     SPACING options to control the behavior of the multi-peak picker.
%
%     The ABSXC option gives the choice to pick positive peaks in the 
%     absolute value of the correlogram or in the unaltered correlogram.  
%     The absolute value of the correlogram can be useful for finding delay
%     times between records/signals that may or may not have opposite 
%     polarity.  A good case for this is when correlating a global set of 
%     body wave phases that were recorded in different portions of the 
%     earthquake's focal sphere.  The unaltered case is useful for times 
%     when there is no polarity difference, for example when correlating 
%     teleseismic phases in a local array.

%     Version History:
%        May   7, 2008 - added convolve stuff? (which didn't work, yay)
%        June 27, 2009 - minor doc fix
%        Sep.  9, 2009 - cleaned up on documentation
%        Oct.  6, 2009 - dropped use of LOGICAL function
%        Oct.  7, 2009 - replaced preallocation zeros w/ nan
%        Dec.  7, 2009 - avoid divide by zero for autocorrelation==0
%        Jan. 27, 2010 - killed convolve option, forced fft/ifft/sum/max to
%                        work on dimension 1, added a simple example, added
%                        a verbose option (requires PRINT_TIME_LEFT)
%        Mar. 12, 2010 - fix for precision issue in normalized case
%        Feb. 15, 2011 - minor doc update
%        Mar. 24, 2012 - minor doc update
%        May  29, 2012 - pow2pad=0 by default
%        Sep. 20, 2012 - pow2pad=1 by default, pow2pad=1+ required due to
%                        wrapping/overlapping of lags for pow2pad=0-
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 20, 2012 at 10:05 GMT

% todo:
% - vectorized gives different lag result if no peaks left (DO NOT CARE)

% at least one input
if(nargin<1)
    error('seizmo:mcxc:badNumInputs','MCXC requires at least one input'); 
end

% defaults
npeaks=0;           % return entire correlograms
spacing=1;          % peaks separated by at least 1 interval
lags=[];            % all possible lags
adjacent=0;         % no adjacent points
normxc=true;        % normalize correlograms
absxc=true;         % take absolute value before picking
vectorxc=false;     % use looped subfunctions
pow2pad=1;          % next power of 2
verbose=true;       % verbosity of messages

% parse and check x and y
if(~isnumeric(x))
    error('seizmo:mcxc:badInput','X must be a numeric array');
end
if(nargin>1 && isnumeric(varargin{1}))
    y=varargin{1};
    varargin(1)=[];
    selfxc=false;
else
    selfxc=true;
end

% parse options
nvarargin=numel(varargin);
if(mod(nvarargin,2))
    error('seizmo:mcxc:unpairedOption',...
        'All options must be ''field''/value pairs');
end
for i=1:2:nvarargin
    % check field is char
    if(~ischar(varargin{i}))
        error('seizmo:mcxc:fieldBad','Option field must be a string');
    end
    
    value=varargin{i+1};
    switch lower(varargin{i})
        case 'npeaks'
            if(~isnumeric(value) || ~isscalar(value) ...
                    || value<0 || (fix(value)-value)~=0)
                error('seizmo:mcxc:badInput',...
                    'NPEAKS must be an integer scalar >= 0'); 
            end
            npeaks=value;
        case 'spacing'
            if(~isnumeric(value) || ~isscalar(value) || value<0)
                error('seizmo:mcxc:badInput',...
                    'SPACING must be an numeric scalar >= 0'); 
            end
            spacing=value;
        case 'lags'
            if(~isnumeric(value) || numel(lags)>2 ...
                    || any((fix(lags)-lags)~=0))
                error('seizmo:mcxc:badInput',['LAGS must be a vector '...
                    'of integers with 0 to 2 elements']);
            end
            lags=value;
        case 'adjacent'
            if(~isnumeric(value) || ~isscalar(value) || value<0)
                error('seizmo:mcxc:badInput',...
                    'ADJACENT must be an numeric scalar >= 0'); 
            end
            adjacent=floor(value);
        case 'normxc'
            if((~isnumeric(value) && ~islogical(value)) ...
                    || ~isscalar(value))
                error('seizmo:mcxc:badInput',...
                    'NORMXC must be a logical scalar!');
            end
            normxc=value;
        case 'absxc'
            if((~isnumeric(value) && ~islogical(value)) ...
                    || ~isscalar(value))
                error('seizmo:mcxc:badInput',...
                    'ABSXC must be a logical scalar!');
            end
            absxc=value;
        case 'vectorxc'
            if((~isnumeric(value) && ~islogical(value)) ...
                    || ~isscalar(value))
                error('seizmo:mcxc:badInput',...
                    'VECTORXC must be a logical scalar!');
            end
            vectorxc=value;
        case 'pow2pad'
            if(~isnumeric(value) || ~isscalar(value) ...
                    || (fix(value)-value)~=0 || value<1)
                error('seizmo:mcxc:badInput',...
                    'POW2PAD must be a positive integer!'); 
            end
            pow2pad=value;
        case 'verbose'
            if((~isnumeric(value) && ~islogical(value)) ...
                    || ~isscalar(value))
                error('seizmo:mcxc:badInput',...
                    'VERBOSE must be a logical scalar!');
            end
            verbose=value;
        otherwise
            error('seizmo:mcxc:unknownOption',...
                'Unknown option: %s',varargin{i});
    end
end

% changing spacing to match subfunction usage
% 0 -> -1 (zero out nothing)
% 0.1 -> 0 (zero out only peak)
% 1 -> 0 (zero out only peak)
% 2.1 -> 2 (zero out peak and 4 neighboring points)
spacing=ceil(spacing-1);

% force lags to be in order
lags=sort(lags);

% detail message
if(verbose)
    disp('Correlating Record(s)')
end

% x or x+y
if(selfxc)
    [cg,lg,pg]=smcxc(x,npeaks,adjacent,spacing,...
        lags,normxc,absxc,vectorxc,pow2pad,verbose);
else
    % size'em up
    sx=size(x); sy=size(y);
    
    % lags
    if(isempty(lags)); lags=[-sx(1)+1 sy(1)-1];
    elseif(isscalar(lags)); lags=sort([-lags lags]);
    end
    
    % length with zeros padding
    nFFT=2^(nextpow2(max([sx(1) sy(1)]))+pow2pad);
    
    % getting transforms and conjugates
    Fx=fft(x,nFFT,1);
    
    % conjugates and autocorrelations
    CFy=conj(fft(y,nFFT,1));
    if(normxc)
        ZLACx=sqrt(sum(x.^2,1));
        ZLACy=sqrt(sum(y.^2,1));

        % avoid divide by zero
        ZLACx(ZLACx==0)=eps;
        ZLACy(ZLACy==0)=eps;
    end
    
    % get peaks
    if(npeaks)
        % vectorize
        if(vectorxc)
            % normalize
            if(normxc)
                % absolute peaks
                if(absxc)
                    [cg,lg,pg]=mcxcmncapv(Fx,CFy,ZLACx,ZLACy,nFFT,sx,sy,npeaks,adjacent,spacing,lags,verbose);
                % positive peaks
                else
                    [cg,lg]=mcxcmncpv(Fx,CFy,ZLACx,ZLACy,nFFT,sx,sy,npeaks,adjacent,spacing,lags,verbose);
                    if(nargout>2); pg=ones(sx(2),sy(2),npeaks,1+2*adjacent); end
                end
            % dont normalize
            else
                % absolute peaks
                if(absxc)
                    [cg,lg,pg]=mcxcmcapv(Fx,CFy,nFFT,sx,sy,npeaks,adjacent,spacing,lags,verbose);
                % positive peaks
                else
                    [cg,lg]=mcxcmcpv(Fx,CFy,nFFT,sx,sy,npeaks,adjacent,spacing,lags,verbose);
                    if(nargout>2); pg=ones(sx(2),sy(2),npeaks,1+2*adjacent); end
                end
            end
        % loop
        else
            % normalize
            if(normxc)
                % absolute peaks
                if(absxc)
                    [cg,lg,pg]=mcxcmncap(Fx,CFy,ZLACx,ZLACy,nFFT,sx,sy,npeaks,adjacent,spacing,lags,verbose);
                % positive peaks
                else
                    [cg,lg]=mcxcmncp(Fx,CFy,ZLACx,ZLACy,nFFT,sx,sy,npeaks,adjacent,spacing,lags,verbose);
                    if(nargout>2); pg=ones(sx(2),sy(2),npeaks,1+2*adjacent); end
                end
            % dont normalize
            else
                % absolute peaks
                if(absxc)
                    [cg,lg,pg]=mcxcmcap(Fx,CFy,nFFT,sx,sy,npeaks,adjacent,spacing,lags,verbose);
                % positive peaks
                else
                    [cg,lg]=mcxcmcp(Fx,CFy,nFFT,sx,sy,npeaks,adjacent,spacing,lags,verbose);
                    if(nargout>2); pg=ones(sx(2),sy(2),npeaks,1+2*adjacent); end
                end
            end
        end
    % get correlograms
    else
        pg=[];
        % vectorize
        if(vectorxc)
            % normalize
            if(normxc)
                [cg,lg]=mcxcncv(Fx,CFy,ZLACx,ZLACy,nFFT,sx,sy,lags,verbose);
            % dont normalize
            else
                [cg,lg]=mcxccv(Fx,CFy,nFFT,sx,sy,lags,verbose);
            end
        % loop
        else
            % normalize
            if(normxc)
                [cg,lg]=mcxcnc(Fx,CFy,ZLACx,ZLACy,nFFT,sx,sy,lags,verbose);
            % dont normalize
            else
                [cg,lg]=mcxcc(Fx,CFy,nFFT,sx,sy,lags,verbose);
            end
        end
    end
end

% Due to autocorrelations being calced in the time domain while cross
% correlations are found in the frequency domain, normalization might be
% a hair off.  This results in values just above/below 1/-1 and they need
% to be shifted to 1/-1.
if(normxc); cg(abs(cg)>1)=sign(cg(abs(cg)>1)); end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                  CORRELOGRAM ONLY MCXC SECTION                  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cg,lags,sh,len]=mcxcc_alloc(nFFT,sx,sy,lagrng)
% allocation for mcxc - looped

% shifts and lags
sh=[(sx(1):-1:1)  (nFFT:-1:nFFT-sy(1)+2)];
lags=(-sx(1)+1:sy(1)-1).';
glags=(lags>=lagrng(1) & lags<=lagrng(2));
lags=lags(glags);
sh=sh(glags);
len=numel(sh);

% preallocate output
cg=nan(sx(2),sy(2),len);

end

function [cg,lags]=mcxcnc(F,CF,ZLACx,ZLACy,nFFT,sx,sy,lagrng,verbose)
%MCXCNC    mcxc 4 normalized correlograms - looped

% detail message
if(verbose)
    nrecs=sx(2)*sy(2);
    print_time_left(0,nrecs);
end

% preallocate just about everything
[cg,lags,shift]=mcxcc_alloc(nFFT,sx,sy,lagrng);

% get correlograms
for i=1:sx(2)
    for j=1:sy(2)
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),[],1,'symmetric');
        cg(i,j,:)=XCF(shift)./(ZLACx(i)*ZLACy(j));
        
        % detail message
        if(verbose)
            print_time_left((i-1)*sy(2)+j,nrecs);
        end
    end
end

end

function [cg,lags]=mcxcc(F,CF,nFFT,sx,sy,lagrng,verbose)
%MCXCC    mcxc 4 correlograms - looped

% detail message
if(verbose)
    nrecs=sx(2)*sy(2);
    print_time_left(0,nrecs);
end

% preallocate just about everything
[cg,lags,shift]=mcxcc_alloc(nFFT,sx,sy,lagrng);

% get correlograms
for i=1:sx(2)
    for j=1:sy(2)
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),[],1,'symmetric');
        cg(i,j,:)=XCF(shift);
        
        % detail message
        if(verbose)
            print_time_left((i-1)*sy(2)+j,nrecs);
        end
    end
end

end

function [cg,lags]=mcxcncv(F,CF,ZLACx,ZLACy,nFFT,sx,sy,lagrng,verbose)
%MCXCNCV    mcxc 4 normalized correlograms - vectorized

% detail message
if(verbose)
    nrecs=sx(2)*sy(2);
    print_time_left(0,nrecs);
end

% allocate
[cg,lags,shift,len]=mcxcc_alloc(nFFT,sx,sy,lagrng);
ZLACy=ZLACy(ones(len,1),:);
cg=permute(cg,[1 3 2]);

% get correlograms
for i=1:sx(2)
    % do all possible pairings against record i
    XCF=ifft(F(:,i*ones(1,sy(2))).*CF,[],1,'symmetric');
    
    % put zero lag at center and normalize
    cg(i,:,:)=XCF(shift,:)./(ZLACx(i)*ZLACy);
        
    % detail message
    if(verbose)
        print_time_left(i*sy(2),nrecs);
    end
end
cg=permute(cg,[1 3 2]);

end

function [cg,lags]=mcxccv(F,CF,nFFT,sx,sy,lagrng,verbose)
%MCXCCV    mcxc 4 normalized correlograms - vectorized

% detail message
if(verbose)
    nrecs=sx(2)*sy(2);
    print_time_left(0,nrecs);
end

% allocate
[cg,lags,shift]=mcxcc_alloc(nFFT,sx,sy,lagrng);
cg=permute(cg,[1 3 2]);

% get correlograms
for i=1:sx(2)
    % do all possible pairings against record i
    XCF=ifft(F(:,i*ones(1,sy(2))).*CF,[],1,'symmetric');
    
    % put zero lag at center
    cg(i,:,:)=XCF(shift,:);
        
    % detail message
    if(verbose)
        print_time_left(i*sy(2),nrecs);
    end
end
cg=permute(cg,[1 3 2]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                      FOR LOOP MCXC SECTION                      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cg,lg,lags,sh,len,range,peak,ptsadj,idxadj]=mcxc_alloc(nFFT,sx,sy,n,adjacent,spacing,lagrng)
% allocation for mcxc - looped

% shifts and lags
sh=[(sx(1):-1:1)  (nFFT:-1:nFFT-sy(1)+2)];
lags=(-sx(1)+1:sy(1)-1).';
glags=(lags>=lagrng(1) & lags<=lagrng(2));
lags=lags(glags);
sh=sh(glags);
len=numel(sh);
range=(-spacing:spacing).';

% preallocate output
cg=nan(sx(2),sy(2),n,1+2*adjacent); lg=cg;

% adjacent
peak=1+adjacent;
ptsadj=[-adjacent:-1 1:adjacent];
idxadj=ptsadj+peak;

end

function [cg,lg]=mcxcmncp(F,CF,ZLACx,ZLACy,nFFT,sx,sy,n,adjacent,spacing,lagrng,verbose)
%MCXCMNCP    mcxc 4 multiple normalized correlogram peaks - looped

% detail message
if(verbose)
    nrecs=sx(2)*sy(2);
    print_time_left(0,nrecs);
end

% preallocate just about everything
[cg,lg,lags,shift,len,range,peak,ptsadj,idxadj]=mcxc_alloc(nFFT,sx,sy,n,adjacent,spacing,lagrng);

% get correlation peaks
for i=1:sx(2)
    for j=1:sy(2)
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),[],1,'symmetric');
        XCFT=XCF(shift)./(ZLACx(i)*ZLACy(j));
        
        % pick peak
        [cg(i,j,1,peak),imxc]=max(XCFT);
        lg(i,j,1,peak)=lags(imxc);
        
        % add adjacent points
        iadj=imxc+ptsadj;
        good=(iadj>0 & iadj<=len);
        cg(i,j,1,idxadj(good))=XCFT(iadj(good));
        lg(i,j,1,idxadj(good))=lags(iadj(good));
        
        % multi-peak picker
        for k=2:n
            % zero previous peak
            izero=imxc+range;
            XCFT(izero(izero>0 & izero<=len))=nan;
            
            % pick peak
            [cg(i,j,k,peak),imxc]=max(XCFT);
            if(isnan(cg(i,j,k,peak))); continue; end
            lg(i,j,k,peak)=lags(imxc);
            
            % add adjacent points
            iadj=imxc+ptsadj;
            good=(iadj>0 & iadj<=len);
            cg(i,j,k,idxadj(good))=XCFT(iadj(good));
            lg(i,j,k,idxadj(good))=lags(iadj(good));
        end
        
        % detail message
        if(verbose)
            print_time_left((i-1)*sy(2)+j,nrecs);
        end
    end
end

end

function [cg,lg]=mcxcmcp(F,CF,nFFT,sx,sy,n,adjacent,spacing,lagrng,verbose)
%MCXCMCP    mcxc 4 multiple correlogram peaks - looped

% detail message
if(verbose)
    nrecs=sx(2)*sy(2);
    print_time_left(0,nrecs);
end

% preallocate just about everything
[cg,lg,lags,shift,len,range,peak,ptsadj,idxadj]=mcxc_alloc(nFFT,sx,sy,n,adjacent,spacing,lagrng);

% get correlation peaks
for i=1:sx(2)
    for j=1:sy(2)
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),[],1,'symmetric');
        XCFT=XCF(shift);
        
        % pick peak
        [cg(i,j,1,peak),imxc]=max(XCFT);
        lg(i,j,1,peak)=lags(imxc);
        
        % add adjacent points
        iadj=imxc+ptsadj;
        good=(iadj>0 & iadj<=len);
        cg(i,j,1,idxadj(good))=XCFT(iadj(good));
        lg(i,j,1,idxadj(good))=lags(iadj(good));
        
        % multi-peak picker
        for k=2:n
            % zero previous peak
            izero=imxc+range;
            XCFT(izero(izero>0 & izero<=len))=nan;
            
            % pick peak
            [cg(i,j,k,peak),imxc]=max(XCFT);
            if(isnan(cg(i,j,k,peak))); continue; end
            lg(i,j,k,peak)=lags(imxc);
            
            % add adjacent points
            iadj=imxc+ptsadj;
            good=(iadj>0 & iadj<=len);
            cg(i,j,k,idxadj(good))=XCFT(iadj(good));
            lg(i,j,k,idxadj(good))=lags(iadj(good));
        end
        
        % detail message
        if(verbose)
            print_time_left((i-1)*sy(2)+j,nrecs);
        end
    end
end

end

function [cg,lg,pg]=mcxcmncap(F,CF,ZLACx,ZLACy,nFFT,sx,sy,n,adjacent,spacing,lagrng,verbose)
%MCXCMNCAP    mcxc 4 multiple normalized correlogram absolute peaks - looped

% detail message
if(verbose)
    nrecs=sx(2)*sy(2);
    print_time_left(0,nrecs);
end

% preallocate just about everything
[cg,lg,lags,shift,len,range,peak,ptsadj,idxadj]=mcxc_alloc(nFFT,sx,sy,n,adjacent,spacing,lagrng);
pg=cg;

% get correlation peaks
for i=1:sx(2)
    for j=1:sy(2)
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),[],1,'symmetric');
        XCFT=XCF(shift)./(ZLACx(i)*ZLACy(j));
        AXCFT=abs(XCFT);
        
        % pick peak
        [cg(i,j,1,peak),imxc]=max(AXCFT);
        lg(i,j,1,peak)=lags(imxc);
        pg(i,j,1,peak)=squeeze(cg(i,j,1,peak))/XCFT(imxc);
        
        % add adjacent points
        iadj=imxc+ptsadj;
        good=(iadj>0 & iadj<=len);
        cg(i,j,1,idxadj(good))=AXCFT(iadj(good));
        lg(i,j,1,idxadj(good))=lags(iadj(good));
        pg(i,j,1,idxadj(good))=squeeze(cg(i,j,1,idxadj(good)))./XCFT(iadj(good));
        
        % multi-peak picker
        for k=2:n
            % zero previous peak
            izero=imxc+range;
            AXCFT(izero(izero>0 & izero<=len))=nan;
            
            % pick peak
            [cg(i,j,k,peak),imxc]=max(AXCFT);
            if(isnan(cg(i,j,k,peak))); continue; end
            lg(i,j,k,peak)=lags(imxc);
            pg(i,j,k,peak)=squeeze(cg(i,j,k,peak))/XCFT(imxc);
            
            % add adjacent points
            iadj=imxc+ptsadj;
            good=(iadj>0 & iadj<=len);
            cg(i,j,k,idxadj(good))=AXCFT(iadj(good));
            lg(i,j,k,idxadj(good))=lags(iadj(good));
            pg(i,j,k,idxadj(good))=squeeze(cg(i,j,k,idxadj(good)))./XCFT(iadj(good));
        end
        
        % detail message
        if(verbose)
            print_time_left((i-1)*sy(2)+j,nrecs);
        end
    end
end

end

function [cg,lg,pg]=mcxcmcap(F,CF,nFFT,sx,sy,n,adjacent,spacing,lagrng,verbose)
%MCXCMCAP    mcxc 4 multiple correlogram absolute peaks - looped

% detail message
if(verbose)
    nrecs=sx(2)*sy(2);
    print_time_left(0,nrecs);
end

% preallocate just about everything
[cg,lg,lags,shift,len,range,peak,ptsadj,idxadj]=mcxc_alloc(nFFT,sx,sy,n,adjacent,spacing,lagrng);
pg=cg;

% get correlation peaks
for i=1:sx(2)
    for j=1:sy(2)
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),[],1,'symmetric');
        XCFT=XCF(shift);
        AXCFT=abs(XCFT);
        
        % pick peak
        [cg(i,j,1,peak),imxc]=max(AXCFT);
        lg(i,j,1,peak)=lags(imxc);
        pg(i,j,1,peak)=squeeze(cg(i,j,1,peak))/XCFT(imxc);
        
        % add adjacent points
        iadj=imxc+ptsadj;
        good=(iadj>0 & iadj<=len);
        cg(i,j,1,idxadj(good))=AXCFT(iadj(good));
        lg(i,j,1,idxadj(good))=lags(iadj(good));
        pg(i,j,1,idxadj(good))=squeeze(cg(i,j,1,idxadj(good)))./XCFT(iadj(good));
        
        % multi-peak picker
        for k=2:n
            % zero previous peak
            izero=imxc+range;
            AXCFT(izero(izero>0 & izero<=len))=nan;
            
            % pick peak
            [cg(i,j,k,peak),imxc]=max(AXCFT);
            if(isnan(cg(i,j,k,peak))); continue; end
            lg(i,j,k,peak)=lags(imxc);
            pg(i,j,k,peak)=squeeze(cg(i,j,k,peak))/XCFT(imxc);
            
            % add adjacent points
            iadj=imxc+ptsadj;
            good=(iadj>0 & iadj<=len);
            cg(i,j,k,idxadj(good))=AXCFT(iadj(good));
            lg(i,j,k,idxadj(good))=lags(iadj(good));
            pg(i,j,k,idxadj(good))=squeeze(cg(i,j,k,idxadj(good)))./XCFT(iadj(good));
        end
        
        % detail message
        if(verbose)
            print_time_left((i-1)*sy(2)+j,nrecs);
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                      VECTORIZED MCXC SECTION                    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cg,lg,lags,sh,len,range,per21,cols,...
    peak,ptsadj,idxadj,nadj,adjcols,cgtemp,lgtemp]=...
    mcxcv_alloc(nFFT,sx,sy,n,adjacent,spacing,lagrng)
% allocation for mcxc - vectorized

% shifts and lags
sh=[(sx(1):-1:1)  (nFFT:-1:nFFT-sy(1)+2)];
lags=(-sx(1)+1:sy(1)-1).';
glags=(lags>=lagrng(1) & lags<=lagrng(2));
lags=lags(glags);
sh=sh(glags);
len=numel(sh);
range=(-spacing:spacing).';
per21=2*spacing+1;

% preallocate output
cg=nan(sx(2),sy(2),n,1+2*adjacent); lg=cg;

% preallocate zeroer
cols=repmat(1:sy(2),per21,1);

% adjacent
peak=1+adjacent;
ptsadj=[-adjacent:-1 1:adjacent].';
idxadj=ptsadj+peak;
nadj=2*adjacent;
adjcols=repmat(1:sy(2),nadj,1);
cgtemp=nan(nadj,sy(2)); 
lgtemp=cgtemp;

end

function [cg,lg]=mcxcmncpv(F,CF,ZLACx,ZLACy,nFFT,sx,sy,n,adjacent,spacing,lagrng,verbose)
%MCXCMNCPV    mcxc 4 multiple normalized correlogram peaks - vectorized

% detail message
if(verbose)
    nrecs=sx(2)*sy(2);
    print_time_left(0,nrecs);
end

% allocate
[cg,lg,lags,shift,len,range,per21,cols,...
    peak,ptsadj,idxadj,nadj,adjcols,cgtemp,lgtemp]=...
    mcxcv_alloc(nFFT,sx,sy,n,adjacent,spacing,lagrng);
ZLACy=ZLACy(ones(len,1),:);

% get correlation peaks
for i=1:sx(2)
    % do all possible pairings against record i
    XCF=ifft(F(:,i*ones(1,sy(2))).*CF,[],1,'symmetric');
    
    % put zero lag at center and normalize
    XCFT=XCF(shift,:)./(ZLACx(i)*ZLACy);
    
    % multi-peak picker
    for j=1:n
        % find peak(s)
        [cg(i,:,j,peak),imxc]=max(XCFT,[],1);
        lg(i,:,j,peak)=lags(imxc);
        
        % add adjacent point(s)
        iadj=imxc(ones(nadj,1),:)+ptsadj(:,ones(1,sy(2)));
        good=(iadj>0 & iadj<=len);
        cgtemp(:,:)=nan;
        lgtemp(:,:)=nan;
        cgtemp(good)=XCFT(sub2ind([len sy(2)],iadj(good),adjcols(good)));
        lgtemp(good)=lags(iadj(good));
        cg(i,:,j,idxadj)=cgtemp.';
        lg(i,:,j,idxadj)=lgtemp.';
        
        % zero out peak(s)
        rows=imxc(ones(per21,1),:)+range(:,ones(1,sy(2)));
        good2=(rows>0 & rows<=len);
        XCFT(sub2ind([len sy(2)],rows(good2),cols(good2)))=nan;
    end
        
    % detail message
    if(verbose)
        print_time_left(i*sy(2),nrecs);
    end
end

end

function [cg,lg]=mcxcmcpv(F,CF,nFFT,sx,sy,n,adjacent,spacing,lagrng,verbose)
%MCXCMCPV    mcxc 4 multiple correlogram peaks - vectorized

% detail message
if(verbose)
    nrecs=sx(2)*sy(2);
    print_time_left(0,nrecs);
end

% allocate
[cg,lg,lags,shift,len,range,per21,cols,...
    peak,ptsadj,idxadj,nadj,adjcols,cgtemp,lgtemp]=...
    mcxcv_alloc(nFFT,sx,sy,n,adjacent,spacing,lagrng);

% get correlation peaks
for i=1:sx(2)
    % do all possible pairings against record i
    XCF=ifft(F(:,i*ones(1,sy(2))).*CF,[],1,'symmetric');
    
    % put zero lag at center
    XCFT=XCF(shift,:);
    
    % multi-peak picker
    for j=1:n
        % find peak(s)
        [cg(i,:,j,peak),imxc]=max(XCFT,[],1);
        lg(i,:,j,peak)=lags(imxc);
        
        % add adjacent point(s)
        iadj=imxc(ones(nadj,1),:)+ptsadj(:,ones(1,sy(2)));
        good=(iadj>0 & iadj<=len);
        cgtemp(:,:)=nan;
        lgtemp(:,:)=nan;
        cgtemp(good)=XCFT(sub2ind([len sy(2)],iadj(good),adjcols(good)));
        lgtemp(good)=lags(iadj(good));
        cg(i,:,j,idxadj)=cgtemp.';
        lg(i,:,j,idxadj)=lgtemp.';
        
        % zero out peak(s)
        rows=imxc(ones(per21,1),:)+range(:,ones(1,sy(2)));
        good2=(rows>0 & rows<=len);
        XCFT(sub2ind([len sy(2)],rows(good2),cols(good2)))=nan;
    end
        
    % detail message
    if(verbose)
        print_time_left(i*sy(2),nrecs);
    end
end

end

function [cg,lg,pg]=mcxcmncapv(F,CF,ZLACx,ZLACy,nFFT,sx,sy,n,adjacent,spacing,lagrng,verbose)
%MCXCMNCAPV    mcxc 4 multiple normalized correlogram absolute peaks - vectorized

% detail message
if(verbose)
    nrecs=sx(2)*sy(2);
    print_time_left(0,nrecs);
end

% allocate
[cg,lg,lags,shift,len,range,per21,cols,...
    peak,ptsadj,idxadj,nadj,adjcols,cgtemp,lgtemp]=...
    mcxcv_alloc(nFFT,sx,sy,n,adjacent,spacing,lagrng);
ZLACy=ZLACy(ones(len,1),:);
pg=cg; pgtemp=cgtemp;

% get correlation peaks
for i=1:sx(2)
    % do all possible pairings against record i
    XCF=ifft(F(:,i*ones(1,sy(2))).*CF,[],1,'symmetric');
    
    % put zero lag at center and make real
    XCFT=XCF(shift,:)./(ZLACx(i)*ZLACy);
    
    % also find absolute values
    AXCFT=abs(XCFT);
    
    % multi-peak picker
    for j=1:n
        % find peak(s)
        [cg(i,:,j,peak),imxc]=max(AXCFT,[],1);
        lg(i,:,j,peak)=lags(imxc);
        pg(i,:,j,peak)=cg(i,:,j,peak)./XCFT((0:sy(2)-1)*len+imxc);
        
        % add adjacent point(s)
        iadj=imxc(ones(nadj,1),:)+ptsadj(:,ones(1,sy(2)));
        good=(iadj>0 & iadj<=len);
        cgtemp(:,:)=nan;
        lgtemp(:,:)=nan;
        pgtemp(:,:)=nan;
        cgtemp(good)=AXCFT(sub2ind([len sy(2)],iadj(good),adjcols(good)));
        lgtemp(good)=lags(iadj(good));
        pgtemp(good)=cgtemp(good)./XCFT(sub2ind([len sy(2)],iadj(good),adjcols(good)));
        cg(i,:,j,idxadj)=cgtemp.';
        lg(i,:,j,idxadj)=lgtemp.';
        pg(i,:,j,idxadj)=pgtemp.';
        
        % zero out peak(s)
        rows=imxc(ones(per21,1),:)+range(:,ones(1,sy(2)));
        good2=(rows>0 & rows<=len);
        AXCFT(sub2ind([len sy(2)],rows(good2),cols(good2)))=nan;
    end
        
    % detail message
    if(verbose)
        print_time_left(i*sy(2),nrecs);
    end
end

end

function [cg,lg,pg]=mcxcmcapv(F,CF,nFFT,sx,sy,n,adjacent,spacing,lagrng,verbose)
%MCXCMCAPV    mcxc 4 multiple correlogram absolute peaks - vectorized

% detail message
if(verbose)
    nrecs=sx(2)*sy(2);
    print_time_left(0,nrecs);
end

% allocate
[cg,lg,lags,shift,len,range,per21,cols,...
    peak,ptsadj,idxadj,nadj,adjcols,cgtemp,lgtemp]=...
    mcxcv_alloc(nFFT,sx,sy,n,adjacent,spacing,lagrng);
pg=cg; pgtemp=cgtemp;

% get correlation peaks
for i=1:sx(2)
    % do all possible pairings against record i
    XCF=ifft(F(:,i*ones(1,sy(2))).*CF,[],1,'symmetric');
    
    % put zero lag at center and make real
    XCFT=XCF(shift,:);
    
    % also find absolute values
    AXCFT=abs(XCFT);
    
    % multi-peak picker
    for j=1:n
        % find peak(s)
        [cg(i,:,j,peak),imxc]=max(AXCFT,[],1);
        lg(i,:,j,peak)=lags(imxc);
        pg(i,:,j,peak)=cg(i,:,j,peak)./XCFT((0:sy(2)-1)*len+imxc);
        
        % add adjacent point(s)
        iadj=imxc(ones(nadj,1),:)+ptsadj(:,ones(1,sy(2)));
        good=(iadj>0 & iadj<=len);
        cgtemp(:,:)=nan;
        lgtemp(:,:)=nan;
        pgtemp(:,:)=nan;
        cgtemp(good)=AXCFT(sub2ind([len sy(2)],iadj(good),adjcols(good)));
        lgtemp(good)=lags(iadj(good));
        pgtemp(good)=cgtemp(good)./XCFT(sub2ind([len sy(2)],iadj(good),adjcols(good)));
        cg(i,:,j,idxadj)=cgtemp.';
        lg(i,:,j,idxadj)=lgtemp.';
        pg(i,:,j,idxadj)=pgtemp.';
        
        % zero out peak(s)
        rows=imxc(ones(per21,1),:)+range(:,ones(1,sy(2)));
        good2=(rows>0 & rows<=len);
        AXCFT(sub2ind([len sy(2)],rows(good2),cols(good2)))=nan;
    end
        
    % detail message
    if(verbose)
        print_time_left(i*sy(2),nrecs);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                          SMCXC SECTION                          %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cv,lv,pv]=smcxc(x,npeaks,adjacent,spacing,lags,normxc,absxc,vectorxc,pow2pad,verbose)
%SMCXC    Self mcxc wrapper

% size'em up
sx=size(x);
txc=(sx(2)^2-sx(2))/2;

% lags
if(isempty(lags)); lags=[-sx(1) sx(1)];
elseif(isscalar(lags)); lags=sort([-lags lags]);
end

% length with zeros padding
nFFT=2^(nextpow2(sx(1))+pow2pad);

% getting transforms and conjugates
F=fft(x,nFFT);

% conjugates and autocorrelations
CF=conj(F);
if(normxc)
    ZLAC=sqrt(sum(x.^2));

    % avoid divide by zero
    ZLAC(ZLAC==0)=eps;
end

% get peaks
if(npeaks)
    % vectorize
    if(vectorxc)
        % normalize
        if(normxc)
            % absolute peaks
            if(absxc)
                [cv,lv,pv]=smcxcmncapv(F,CF,ZLAC,nFFT,sx,npeaks,adjacent,spacing,lags,verbose);
            % positive peaks
            else
                [cv,lv]=smcxcmncpv(F,CF,ZLAC,nFFT,sx,npeaks,adjacent,spacing,lags,verbose);
                if(nargout>2); pv=ones(txc,1,npeaks,1+2*adjacent); end
            end
        % dont normalize
        else
            % absolute peaks
            if(absxc)
                [cv,lv,pv]=smcxcmcapv(F,CF,nFFT,sx,npeaks,adjacent,spacing,lags,verbose);
            % positive peaks
            else
                [cv,lv]=smcxcmcpv(F,CF,nFFT,sx,npeaks,adjacent,spacing,lags,verbose);
                if(nargout>2); pv=ones(txc,1,npeaks,1+2*adjacent); end
            end
        end
    % loop
    else
        % normalize
        if(normxc)
            % absolute peaks
            if(absxc)
                [cv,lv,pv]=smcxcmncap(F,CF,ZLAC,nFFT,sx,npeaks,adjacent,spacing,lags,verbose);
            % positive peaks
            else
                [cv,lv]=smcxcmncp(F,CF,ZLAC,nFFT,sx,npeaks,adjacent,spacing,lags,verbose);
                if(nargout>2); pv=ones(txc,1,npeaks,1+2*adjacent); end
            end
        % dont normalize
        else
            % absolute peaks
            if(absxc)
                [cv,lv,pv]=smcxcmcap(F,CF,nFFT,sx,npeaks,adjacent,spacing,lags,verbose);
            % positive peaks
            else
                [cv,lv]=smcxcmcp(F,CF,nFFT,sx,npeaks,adjacent,spacing,lags,verbose);
                if(nargout>2); pv=ones(txc,1,npeaks,1+2*adjacent); end
            end
        end
    end
% get correlograms
else
    pv=[];
    % vectorize
    if(vectorxc)
        % normalize
        if(normxc)
            [cv,lv]=smcxcncv(F,CF,ZLAC,nFFT,sx,lags,verbose);
        % dont normalize
        else
            [cv,lv]=smcxccv(F,CF,nFFT,sx,lags,verbose);
        end
    % loop
    else
        % normalize
        if(normxc)
            [cv,lv]=smcxcnc(F,CF,ZLAC,nFFT,sx,lags,verbose);
        % dont normalize
        else
            [cv,lv]=smcxcc(F,CF,nFFT,sx,lags,verbose);
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 CORRELOGRAM ONLY SMCXC SECTION                  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cv,lags,sh,c,txc]=smcxcc_alloc(nFFT,sx,lagrng)
% allocation for smcxcc - looped

% shifts and lags
sh=[(sx(1):-1:1)  (nFFT:-1:nFFT-sx(1)+2)];
lags=(-sx(1)+1:sx(1)-1).';
glags=(lags>=lagrng(1) & lags<=lagrng(2));
lags=lags(glags);
sh=sh(glags);
len=numel(sh);
ilen=sx(2)-1:-1:1;
c=[0 cumsum(ilen)];
txc=(sx(2)^2-sx(2))/2;

% preallocate output
cv=nan(txc,1,len);

end

function [cv,lags]=smcxcnc(F,CF,ZLAC,nFFT,sx,lagrng,verbose)
%SMCXCNC    smcxc 4 normalized correlograms - looped

% preallocate just about everything
[cv,lags,shift,c,nrecs]=smcxcc_alloc(nFFT,sx,lagrng);

% detail message
if(verbose)
    print_time_left(0,nrecs);
end

% get correlograms
for i=1:sx(2)-1
    for j=i+1:sx(2)
        k=c(i)+j-i;
        
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),[],1,'symmetric');
        cv(k,1,:)=XCF(shift)./(ZLAC(i)*ZLAC(j));
        
        % detail message
        if(verbose)
            print_time_left(k,nrecs);
        end
    end
end

end

function [cv,lags]=smcxcc(F,CF,nFFT,sx,lagrng,verbose)
%SMCXCC    smcxc 4 correlograms - looped

% preallocate just about everything
[cv,lags,shift,c,nrecs]=smcxcc_alloc(nFFT,sx,lagrng);

% detail message
if(verbose)
    print_time_left(0,nrecs);
end

% get correlograms
for i=1:sx(2)-1
    for j=i+1:sx(2)
        k=c(i)+j-i;
        
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),[],1,'symmetric');
        cv(k,1,:)=XCF(shift);
        
        % detail message
        if(verbose)
            print_time_left(k,nrecs);
        end
    end
end

end

function [cv,XCF,lags,sh,len,ilen,c,c2,txc]=smcxccv_alloc(nFFT,sx,lagrng)
% allocation for smcxcc - vectorized

% preallocate shifts and lags
sh=[(sx(1):-1:1)  (nFFT:-1:nFFT-sx(1)+2)];
lags=(-sx(1)+1:sx(1)-1).';
glags=(lags>=lagrng(1) & lags<=lagrng(2));
lags=lags(glags);
sh=sh(glags);
len=numel(sh);
ilen=sx(2)-1:-1:1;
c2=cumsum(ilen);
c=[1 c2+1];
txc=(sx(2)^2-sx(2))/2;

% preallocate output
cv=nan(txc,1,len);

% preallocate correlograms
XCF=nan(nFFT,sx(2));

end

function [cv,lags]=smcxcncv(F,CF,ZLAC,nFFT,sx,lagrng,verbose)
%SMCXCNCV    smcxc 4 normalized correlograms - vectorized

% preallocate just about everything
[cv,XCF,lags,shift,len,ilen,c,c2,nrecs]=smcxccv_alloc(nFFT,sx,lagrng);
ZLAC2=ZLAC(ones(len,1),:);
cv=permute(cv,[3 2 1]);

% detail message
if(verbose)
    print_time_left(0,nrecs);
end

% get correlograms
for i=1:sx(2)-1
    % indicies
    i2=i+1:sx(2);
    
    % do all possible pairings against record i
    XCF(:,i2)=ifft(F(:,i*ones(1,ilen(i))).*CF(:,i2),[],1,'symmetric');
    
    % put zero lag at center & normalize
    cv(:,1,c(i):c2(i))=XCF(shift,i2)./(ZLAC(i)*ZLAC2(:,i2));
    
    % detail message
    if(verbose)
        print_time_left(c2(i),nrecs);
    end
end
cv=permute(cv,[3 2 1]);

end

function [cv,lags]=smcxccv(F,CF,nFFT,sx,lagrng,verbose)
%SMCXCCV    smcxc 4 correlograms - vectorized

% preallocate just about everything
[cv,XCF,lags,shift,len,ilen,c,c2,nrecs]=smcxccv_alloc(nFFT,sx,lagrng);
cv=permute(cv,[3 2 1]);

% detail message
if(verbose)
    print_time_left(0,nrecs);
end

% get correlograms
for i=1:sx(2)-1
    % indicies
    i2=i+1:sx(2);
    
    % do all possible pairings against record i
    XCF(:,i2)=ifft(F(:,i*ones(1,ilen(i))).*CF(:,i2),[],1,'symmetric');
    
    % put zero lag at center
    cv(:,1,c(i):c2(i))=XCF(shift,i2);
    
    % detail message
    if(verbose)
        print_time_left(c2(i),nrecs);
    end
end
cv=permute(cv,[3 2 1]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                       LOOPED SMCXC SECTION                      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cv,lv,lags,sh,len,range,per21,c,peak,ptsadj,idxadj,txc]=smcxc_alloc(nFFT,sx,n,adjacent,spacing,lagrng)
% allocation for smcxc - looped

% preallocate shifts and lags
sh=[(sx(1):-1:1)  (nFFT:-1:nFFT-sx(1)+2)];
lags=(-sx(1)+1:sx(1)-1).';
glags=(lags>=lagrng(1) & lags<=lagrng(2));
lags=lags(glags);
sh=sh(glags);
len=numel(sh);
range=(-spacing:spacing).';
per21=2*spacing+1;
ilen=sx(2)-1:-1:1;
c=[0 cumsum(ilen)];
txc=(sx(2)^2-sx(2))/2;

% preallocate output
cv=nan(txc,1,n,1+2*adjacent); lv=cv;

% adjacent
peak=1+adjacent;
ptsadj=[-adjacent:-1 1:adjacent];
idxadj=ptsadj+peak;

end

function [cv,lv]=smcxcmncp(F,CF,ZLAC,nFFT,sx,n,adjacent,spacing,lagrng,verbose)
%SMCXCMNCP    smcxc 4 multiple normalized correlogram peaks - looped

% preallocate just about everything
[cv,lv,lags,shift,len,range,per21,c,peak,ptsadj,idxadj,nrecs]=...
    smcxc_alloc(nFFT,sx,n,adjacent,spacing,lagrng);

% detail message
if(verbose)
    print_time_left(0,nrecs);
end

% get correlation peaks
for i=1:sx(2)-1
    for j=i+1:sx(2)
        k=c(i)+j-i;
        
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),[],1,'symmetric');
        XCFT=XCF(shift)./(ZLAC(i)*ZLAC(j));
        
        % pick peak
        [cv(k,1,1,peak),imxc]=max(XCFT);
        lv(k,1,1,peak)=lags(imxc);
        
        % add adjacent points
        iadj=imxc+ptsadj;
        good=(iadj>0 & iadj<=len);
        cv(k,1,1,idxadj(good))=XCFT(iadj(good));
        lv(k,1,1,idxadj(good))=lags(iadj(good));
        
        % multi-peak picker
        for p=2:n
            % zero previous peak
            izero=imxc+range;
            XCFT(izero(izero>0 & izero<=len))=nan;
            
            % pick peak
            [cv(k,1,p,peak),imxc]=max(XCFT);
            if(isnan(cv(k,1,p,peak))); continue; end
            lv(k,1,p,peak)=lags(imxc);
            
            % add adjacent points
            iadj=imxc+ptsadj;
            good=(iadj>0 & iadj<=len);
            cv(k,1,p,idxadj(good))=XCFT(iadj(good));
            lv(k,1,p,idxadj(good))=lags(iadj(good));
        end
        
        % detail message
        if(verbose)
            print_time_left(k,nrecs);
        end
    end
end

end

function [cv,lv]=smcxcmcp(F,CF,nFFT,sx,n,adjacent,spacing,lagrng,verbose)
%SMCXCMCP    smcxc 4 multiple correlogram peaks - looped

% preallocate just about everything
[cv,lv,lags,shift,len,range,per21,c,peak,ptsadj,idxadj,nrecs]=...
    smcxc_alloc(nFFT,sx,n,adjacent,spacing,lagrng);

% detail message
if(verbose)
    print_time_left(0,nrecs);
end

% get correlation peaks
for i=1:sx(2)-1
    for j=i+1:sx(2)
        k=c(i)+j-i;
        
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),[],1,'symmetric');
        XCFT=XCF(shift);
        
        % pick peak
        [cv(k,1,1,peak),imxc]=max(XCFT);
        lv(k,1,1,peak)=lags(imxc);
        
        % add adjacent points
        iadj=imxc+ptsadj;
        good=(iadj>0 & iadj<=len);
        cv(k,1,1,idxadj(good))=XCFT(iadj(good));
        lv(k,1,1,idxadj(good))=lags(iadj(good));
        
        % multi-peak picker
        for p=2:n
            % zero previous peak
            izero=imxc+range;
            XCFT(izero(izero>0 & izero<=len))=nan;
            
            % pick peak
            [cv(k,1,p,peak),imxc]=max(XCFT);
            if(isnan(cv(k,1,p,peak))); continue; end
            lv(k,1,p,peak)=lags(imxc);
            
            % add adjacent points
            iadj=imxc+ptsadj;
            good=(iadj>0 & iadj<=len);
            cv(k,1,p,idxadj(good))=XCFT(iadj(good));
            lv(k,1,p,idxadj(good))=lags(iadj(good));
        end
        
        % detail message
        if(verbose)
            print_time_left(k,nrecs);
        end
    end
end

end

function [cv,lv,pv]=smcxcmncap(F,CF,ZLAC,nFFT,sx,n,adjacent,spacing,lagrng,verbose)
%SMCXCMNCAP    smcxc 4 multiple normalized correlogram absolute peaks - looped

% preallocate just about everything
[cv,lv,lags,shift,len,range,per21,c,peak,ptsadj,idxadj,nrecs]=...
    smcxc_alloc(nFFT,sx,n,adjacent,spacing,lagrng);
pv=nan(nrecs,1,n,1+2*adjacent);

% detail message
if(verbose)
    print_time_left(0,nrecs);
end

% get correlation peaks
for i=1:sx(2)-1
    for j=i+1:sx(2)
        k=c(i)+j-i;
        
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),[],1,'symmetric');
        XCFT=XCF(shift)./(ZLAC(i)*ZLAC(j));
        AXCFT=abs(XCFT);
        
        % pick peak
        [cv(k,1,1,peak),imxc]=max(AXCFT);
        lv(k,1,1,peak)=lags(imxc);
        pv(k,1,1,peak)=cv(k,1,1,peak)/XCFT(imxc);
        
        % add adjacent points
        iadj=imxc+ptsadj;
        good=(iadj>0 & iadj<=len);
        cv(k,1,1,idxadj(good))=AXCFT(iadj(good));
        lv(k,1,1,idxadj(good))=lags(iadj(good));
        pv(k,1,1,idxadj(good))=squeeze(cv(k,1,1,idxadj(good)))./XCFT(iadj(good));
        
        % multi-peak picker
        for p=2:n
            % zero previous peak
            izero=imxc+range;
            AXCFT(izero(izero>0 & izero<=len))=nan;
            
            % pick peak
            [cv(k,1,p,peak),imxc]=max(AXCFT);
            if(isnan(cv(k,1,p,peak))); continue; end
            lv(k,1,p,peak)=lags(imxc);
            pv(k,1,p,peak)=cv(k,1,p,peak)/XCFT(imxc);
            
            % add adjacent points
            iadj=imxc+ptsadj;
            good=(iadj>0 & iadj<=len);
            cv(k,1,p,idxadj(good))=AXCFT(iadj(good));
            lv(k,1,p,idxadj(good))=lags(iadj(good));
            pv(k,1,p,idxadj(good))=squeeze(cv(k,1,p,idxadj(good)))./XCFT(iadj(good));
        end
        
        % detail message
        if(verbose)
            print_time_left(k,nrecs);
        end
    end
end

end

function [cv,lv,pv]=smcxcmcap(F,CF,nFFT,sx,n,adjacent,spacing,lagrng,verbose)
%SMCXCMCAP    smcxc 4 multiple correlogram absolute peaks - looped

% preallocate just about everything
[cv,lv,lags,shift,len,range,per21,c,peak,ptsadj,idxadj,nrecs]=...
    smcxc_alloc(nFFT,sx,n,adjacent,spacing,lagrng);
pv=nan(nrecs,1,n,1+2*adjacent);

% detail message
if(verbose)
    print_time_left(0,nrecs);
end

% get correlation peaks
for i=1:sx(2)-1
    for j=i+1:sx(2)
        k=c(i)+j-i;
        
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),[],1,'symmetric');
        XCFT=XCF(shift);
        AXCFT=abs(XCFT);
        
        % pick peak
        [cv(k,1,1,peak),imxc]=max(AXCFT);
        lv(k,1,1,peak)=lags(imxc);
        pv(k,1,1,peak)=cv(k,1,1,peak)/XCFT(imxc);
        
        % add adjacent points
        iadj=imxc+ptsadj;
        good=(iadj>0 & iadj<=len);
        cv(k,1,1,idxadj(good))=AXCFT(iadj(good));
        lv(k,1,1,idxadj(good))=lags(iadj(good));
        pv(k,1,1,idxadj(good))=squeeze(cv(k,1,1,idxadj(good)))./XCFT(iadj(good));
        
        % multi-peak picker
        for p=2:n
            % zero previous peak
            izero=imxc+range;
            AXCFT(izero(izero>0 & izero<=len))=nan;
            
            % pick peak
            [cv(k,1,p,peak),imxc]=max(AXCFT);
            if(isnan(cv(k,1,p,peak))); continue; end
            lv(k,1,p,peak)=lags(imxc);
            pv(k,1,p,peak)=cv(k,1,p,peak)/XCFT(imxc);
            
            % add adjacent points
            iadj=imxc+ptsadj;
            good=(iadj>0 & iadj<=len);
            cv(k,1,p,idxadj(good))=AXCFT(iadj(good));
            lv(k,1,p,idxadj(good))=lags(iadj(good));
            pv(k,1,p,idxadj(good))=squeeze(cv(k,1,p,idxadj(good)))./XCFT(iadj(good));
        end
        
        % detail message
        if(verbose)
            print_time_left(k,nrecs);
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                    VECTORIZED SMCXC SECTION                     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cv,lv,XCF,XCFT,lags,sh,len,range,per21,ilen,c,c2,imxc,rows,cols,...
    good2,peak,ptsadj,idxadj,nadj,adjcols,cgtemp,lgtemp,iadj,good,txc]=...
    smcxcv_alloc(nFFT,sx,n,adjacent,spacing,lagrng)
% allocation for smcxc - vectorized

% preallocate shifts and lags
sh=[(sx(1):-1:1)  (nFFT:-1:nFFT-sx(1)+2)];
lags=(-sx(1)+1:sx(1)-1).';
glags=(lags>=lagrng(1) & lags<=lagrng(2));
lags=lags(glags);
sh=sh(glags);
len=numel(sh);
range=(-spacing:spacing).';
per21=2*spacing+1;
ilen=sx(2)-1:-1:1;
c2=cumsum(ilen);
c=[1 c2+1];
txc=(sx(2)^2-sx(2))/2;
imxc=nan(1,sx(2));

% preallocate output
cv=nan(txc,1,n,1+2*adjacent); lv=cv;

% preallocate correlograms
XCF=nan(nFFT,sx(2));
XCFT=nan(len,sx(2));

% preallocate peak zeroer
rows=nan(per21,sx(2));
cols=repmat(1:sx(2),per21,1);
good2=false(size(rows));

% adjacent
peak=1+adjacent;
ptsadj=[-adjacent:-1 1:adjacent].';
idxadj=ptsadj+peak;
nadj=2*adjacent;
adjcols=repmat(1:sx(2),nadj,1);
cgtemp=nan(nadj,sx(2)); 
lgtemp=cgtemp;
good=false(size(cgtemp));
iadj=cgtemp;

end

function [cv,lv]=smcxcmncpv(F,CF,ZLAC,nFFT,sx,n,adjacent,spacing,lagrng,verbose)
%SMCXCMNCPV    smcxc 4 multiple normalized correlogram peaks - vectorized

% preallocate just about everything
[cv,lv,XCF,XCFT,lags,shift,len,range,per21,ilen,c,c2,imxc,rows,cols,good2,...
    peak,ptsadj,idxadj,nadj,adjcols,cgtemp,lgtemp,iadj,good,nrecs]=...
    smcxcv_alloc(nFFT,sx,n,adjacent,spacing,lagrng);
ZLAC2=ZLAC(ones(len,1),:);

% detail message
if(verbose)
    print_time_left(0,nrecs);
end

% get correlation peaks
for i=1:sx(2)-1
    % indicies
    i2=i+1:sx(2);
    
    % do all possible pairings against record i
    XCF(:,i2)=ifft(F(:,i*ones(1,ilen(i))).*CF(:,i2),[],1,'symmetric');
    
    % put zero lag at center & normalize
    XCFT(:,i2)=XCF(shift,i2)./(ZLAC(i)*ZLAC2(:,i2));
    
    % multi-peak picker
    for j=1:n
        % find peak(s)
        %[cg(i,:,j,peak),imxc]=max(XCFT);
        %lg(i,:,j,peak)=lags(imxc);
        [cv(c(i):c2(i),1,j,peak),imxc(i2)]=max(XCFT(:,i2),[],1);
        lv(c(i):c2(i),1,j,peak)=lags(imxc(i2));
        
        % add adjacent point(s)
        %iadj=imxc(ones(nadj,1),:)+ptsadj(:,ones(1,sy(2)));
        %good=(iadj>0 & iadj<=len);
        %cgtemp(:,:)=nan;
        %lgtemp(:,:)=nan;
        %cgtemp(good)=XCFT(sub2ind([len sy(2)],iadj(good),adjcols(good)));
        %lgtemp(good)=lags(iadj(good));
        %cg(i,:,j,idxadj)=cgtemp.';
        %lg(i,:,j,idxadj)=lgtemp.';
        iadj(:,i2)=imxc(ones(nadj,1),i2)+ptsadj(:,ones(1,ilen(i)));
        good(:,:)=false; cgtemp(:,:)=nan; lgtemp(:,:)=nan;
        good(:,i2)=(iadj(:,i2)>0 & iadj(:,i2)<=len);
        cgtemp(good)=XCFT(sub2ind([len sx(2)],iadj(good),adjcols(good)));
        lgtemp(good)=lags(iadj(good));
        cv(c(i):c2(i),:,j,idxadj)=cgtemp(:,i2).';
        lv(c(i):c2(i),:,j,idxadj)=lgtemp(:,i2).';
        
        % zero out peak(s)
        %rows=imxc(ones(per21,1),:)+range(:,ones(1,sy(2)));
        %good=(rows>0 & rows<=len);
        %XCFT(sub2ind([len sy(2)],rows(good),cols(good)))=nan;
        rows(:,i2)=imxc(ones(per21,1),i2)+range(:,ones(1,ilen(i)));
        good2(:,:)=false; good2(:,i2)=(rows(:,i2)>0 & rows(:,i2)<=len);
        XCFT(sub2ind([len sx(2)],rows(good2),cols(good2)))=nan;
    end
    
    % detail message
    if(verbose)
        print_time_left(c2(i),nrecs);
    end
end

end

function [cv,lv]=smcxcmcpv(F,CF,nFFT,sx,n,adjacent,spacing,lagrng,verbose)
%SMCXCMCPV    smcxc 4 multiple correlogram peaks - vectorized

% preallocate just about everything
[cv,lv,XCF,XCFT,lags,shift,len,range,per21,ilen,c,c2,imxc,rows,cols,good2,...
    peak,ptsadj,idxadj,nadj,adjcols,cgtemp,lgtemp,iadj,good,nrecs]=...
    smcxcv_alloc(nFFT,sx,n,adjacent,spacing,lagrng);

% detail message
if(verbose)
    print_time_left(0,nrecs);
end

% get correlation peaks
for i=1:sx(2)-1
    % indicies
    i2=i+1:sx(2);
    
    % do all possible pairings against record i
    XCF(:,i2)=ifft(F(:,i*ones(1,ilen(i))).*CF(:,i2),[],1,'symmetric');
    
    % put zero lag at center
    XCFT(:,i2)=XCF(shift,i2);
    
    % multi-peak picker
    for j=1:n
        % find peak(s)
        [cv(c(i):c2(i),1,j,peak),imxc(i2)]=max(XCFT(:,i2),[],1);
        lv(c(i):c2(i),1,j,peak)=lags(imxc(i2));
        
        % add adjacent point(s)
        iadj(:,i2)=imxc(ones(nadj,1),i2)+ptsadj(:,ones(1,ilen(i)));
        good(:,:)=false; cgtemp(:,:)=nan; lgtemp(:,:)=nan;
        good(:,i2)=(iadj(:,i2)>0 & iadj(:,i2)<=len);
        cgtemp(good)=XCFT(sub2ind([len sx(2)],iadj(good),adjcols(good)));
        lgtemp(good)=lags(iadj(good));
        cv(c(i):c2(i),:,j,idxadj)=cgtemp(:,i2).';
        lv(c(i):c2(i),:,j,idxadj)=lgtemp(:,i2).';
        
        % zero out peak(s)
        rows(:,i2)=imxc(ones(per21,1),i2)+range(:,ones(1,ilen(i)));
        good2(:,:)=false; good2(:,i2)=(rows(:,i2)>0 & rows(:,i2)<=len);
        XCFT(sub2ind([len sx(2)],rows(good2),cols(good2)))=nan;
    end
    
    % detail message
    if(verbose)
        print_time_left(c2(i),nrecs);
    end
end

end

function [cv,lv,pv]=smcxcmncapv(F,CF,ZLAC,nFFT,sx,n,adjacent,spacing,lagrng,verbose)
%SMCXCMNCAPV    smcxc 4 multiple normalized correlogram absolute peaks - vectorized

% preallocate just about everything
[cv,lv,XCF,XCFT,lags,shift,len,range,per21,ilen,c,c2,imxc,rows,cols,good2,...
    peak,ptsadj,idxadj,nadj,adjcols,cgtemp,lgtemp,iadj,good,nrecs]=...
    smcxcv_alloc(nFFT,sx,n,adjacent,spacing,lagrng);
ZLAC2=ZLAC(ones(len,1),:);
pv=cv; AXCFT=XCFT; pgtemp=cgtemp;

% detail message
if(verbose)
    print_time_left(0,nrecs);
end

% get correlation peaks
for i=1:sx(2)-1
    % indicies
    i2=i+1:sx(2);
    
    % do all possible pairings against record i
    XCF(:,i2)=ifft(F(:,i*ones(1,ilen(i))).*CF(:,i2),[],1,'symmetric');
    
    % put zero lag at center & normalize
    XCFT(:,i2)=XCF(shift,i2)./(ZLAC(i)*ZLAC2(:,i2));
    
    % also find absolute values
    AXCFT(:,i2)=abs(XCFT(:,i2));
    
    % multi-peak picker
    for j=1:n
        % find peak(s)
        [cv(c(i):c2(i),1,j,peak),imxc(i2)]=max(AXCFT(:,i2),[],1);
        lv(c(i):c2(i),1,j,peak)=lags(imxc(i2));
        pv(c(i):c2(i),1,j,peak)=cv(c(i):c2(i),1,j,peak)./XCFT(imxc(i2)+(i2-1)*len).';
        
        % add adjacent point(s)
        iadj(:,i2)=imxc(ones(nadj,1),i2)+ptsadj(:,ones(1,ilen(i)));
        good(:,:)=false; cgtemp(:,:)=nan; lgtemp(:,:)=nan; pgtemp(:,:)=nan;
        good(:,i2)=(iadj(:,i2)>0 & iadj(:,i2)<=len);
        cgtemp(good)=AXCFT(sub2ind([len sx(2)],iadj(good),adjcols(good)));
        lgtemp(good)=lags(iadj(good));
        pgtemp(good)=cgtemp(good)./XCFT(sub2ind([len sx(2)],iadj(good),adjcols(good)));
        cv(c(i):c2(i),:,j,idxadj)=cgtemp(:,i2).';
        lv(c(i):c2(i),:,j,idxadj)=lgtemp(:,i2).';
        pv(c(i):c2(i),:,j,idxadj)=pgtemp(:,i2).';
        
        % zero out peak(s)
        rows(:,i2)=imxc(ones(per21,1),i2)+range(:,ones(1,ilen(i)));
        good2(:,:)=false; good2(:,i2)=(rows(:,i2)>0 & rows(:,i2)<=len);
        AXCFT(sub2ind([len sx(2)],rows(good2),cols(good2)))=nan;
    end
    
    % detail message
    if(verbose)
        print_time_left(c2(i),nrecs);
    end
end

end

function [cv,lv,pv]=smcxcmcapv(F,CF,nFFT,sx,n,adjacent,spacing,lagrng,verbose)
%SMCXCMCAPV    smcxc 4 multiple correlogram absolute peaks - vectorized

% preallocate just about everything
[cv,lv,XCF,XCFT,lags,shift,len,range,per21,ilen,c,c2,imxc,rows,cols,good2,...
    peak,ptsadj,idxadj,nadj,adjcols,cgtemp,lgtemp,iadj,good,nrecs]=...
    smcxcv_alloc(nFFT,sx,n,adjacent,spacing,lagrng);
pv=cv; AXCFT=XCFT; pgtemp=cgtemp;

% detail message
if(verbose)
    print_time_left(0,nrecs);
end

% get correlation peaks
for i=1:sx(2)-1
    % indicies
    i2=i+1:sx(2);
    
    % do all possible pairings against record i
    XCF(:,i2)=ifft(F(:,i*ones(1,ilen(i))).*CF(:,i2),[],1,'symmetric');
    
    % put zero lag at center & normalize
    XCFT(:,i2)=XCF(shift,i2);
    
    % also find absolute values
    AXCFT(:,i2)=abs(XCFT(:,i2));
    
    % multi-peak picker
    for j=1:n
        % find peak(s)
        [cv(c(i):c2(i),1,j,peak),imxc(i2)]=max(AXCFT(:,i2),[],1);
        lv(c(i):c2(i),1,j,peak)=lags(imxc(i2));
        pv(c(i):c2(i),1,j,peak)=cv(c(i):c2(i),1,j,peak)./XCFT(imxc(i2)+(i2-1)*len).';
        
        % add adjacent point(s)
        iadj(:,i2)=imxc(ones(nadj,1),i2)+ptsadj(:,ones(1,ilen(i)));
        good(:,:)=false; cgtemp(:,:)=nan; lgtemp(:,:)=nan; pgtemp(:,:)=nan;
        good(:,i2)=(iadj(:,i2)>0 & iadj(:,i2)<=len);
        cgtemp(good)=AXCFT(sub2ind([len sx(2)],iadj(good),adjcols(good)));
        lgtemp(good)=lags(iadj(good));
        pgtemp(good)=cgtemp(good)./XCFT(sub2ind([len sx(2)],iadj(good),adjcols(good)));
        cv(c(i):c2(i),:,j,idxadj)=cgtemp(:,i2).';
        lv(c(i):c2(i),:,j,idxadj)=lgtemp(:,i2).';
        pv(c(i):c2(i),:,j,idxadj)=pgtemp(:,i2).';
        
        % zero out peak(s)
        rows(:,i2)=imxc(ones(per21,1),i2)+range(:,ones(1,ilen(i)));
        good2(:,:)=false; good2(:,i2)=(rows(:,i2)>0 & rows(:,i2)<=len);
        AXCFT(sub2ind([len sx(2)],rows(good2),cols(good2)))=nan;
    end
    
    % detail message
    if(verbose)
        print_time_left(c2(i),nrecs);
    end
end

end
