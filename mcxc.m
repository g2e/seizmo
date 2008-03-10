function [cg,lg,pg]=mcxc(x,y,npeaks,spacing,lags,normxc,absxc,vectorxc,pow2pad)
%MCXC    Return correlogram peaks from multi-channel X correlation
%
%    Description: Performs multi-channel cross correlation of records
%     distributed in the columns of X against records distributed in the
%     columns of Y and then finds and returns info on the strongest NPEAKS 
%     peaks in each cross correlogram that are spaced at least SPACING 
%     sample intervals apart and are in the lag range LAGS.  LAGS can be a
%     one (symmetric window) or two (non-symmetric window) element vector 
%     indicating the range of lags (in samples) to search. To normalize the 
%     cross correlograms set NORMXC to true.  To pick both positive and 
%     negative peaks in the correlogram set ABSXC to true.  To use the 
%     vectorized cross correlation subfunctions (WHICH ARE TYPICALLY SLOWER
%     due to Matlab's JIT compiler preferring small looped operations over 
%     chunky vectorized operations) set VECTORXC to true.  POW2PAD controls
%     the length padding for spectral operations and so controls the 
%     frequency resolution.  POW2PAD may be any integer. Negative POW2PAD 
%     will truncate X and Y, POW2PAD==0 will minimally pad X and Y, and 
%     positive POW2PAD pad X and Y even more. Using POW2PAD outside of 0-3 
%     is NOT recommended.
%
%     Defaults:                                Range:
%      - NORMXC   == FALSE                      TRUE / FALSE
%      - ABSXC    == FALSE                      TRUE / FALSE
%      - VECTORXC == FALSE                      TRUE / FALSE
%      - NPEAKS   == 1                          0+
%      - SPACING  == 1                          1+
%      - LAGS     == [-size(X,1)+1 size(Y,1)-1] -size(X,1)+1 to size(Y,1)-1
%      - POW2PAD  == 1                          INTEGERS (0-3 RECOMMENDED)
%
%     Empty Y: When Y=[], MCXC cross correlates X against itself.  Output
%      is changed to save space and computation time.
%     
%     NPEAKS==0: When NPEAKS==0, MCXC returns the correlograms rather than
%      peaks.  Note that this can easily generate a massive array that may 
%      exceed your platform's capacity.  For this option ABSXC, NPEAKS, 
%      and SPACING are ignored, PG=[], LG is the vector of integer lags, 
%      and CG has correlograms oriented down the third dimension (CG(1,1,:) 
%      would give the first correlogram with lags given by LG).
%
%     Output:
%
%      CG,LG,PG
%                   SIZE                          Y=[] SIZE
%
%        [size(X,2),size(Y,2),NPEAKS] [size(X,2)*(size(X,2)-1)/2,1,NPEAKS]
%
%      CG: Contains the correlation values at peaks ordered from strongest
%       to weakest so that CG(:,:,1) contains the strongest peak from every
%       correlogram, CG(:,:,2) contains the second strongest peaks, etc.  
%       Rows correspond to the records in X and columns to Y so that
%       CG(4,7,3) would correspond to the 3rd strongest peak in the
%       correlogram from cross correlating X(:,4) with Y(:,7);
%
%      LG: Lags corresponding to peak correlation values given in CG.
%
%      PG: Signs of peaks in CG.  Only useful when ABSXC is TRUE.
%    
%    Usage: [CG,LG,PG]=mcxc(X,Y,NPEAKS,SPACING,NORMXC,ABSXC,...
%                                   VECTORXC,POW2PAD)
%
%    See also: xcorr, fftfilt (Signal Processing Toolbox)

%    Seismology Specific Notes:
%     The idea behind mcxc is that it is a more robust way of determining
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
%     This particular version implements a multiple peak picker, which is 
%     particularly useful for dealing with situations where phase-skipping 
%     is an issue, such as in noisy and/or narrow-band applications.  Use
%     the NPEAKS and SPACING options to control the behavior of the 
%     multi-picker.
%
%     The ABSXC option gives the choice to look for the maximum positive 
%     peaks in the absolute value of the correlogram or in the unaltered 
%     correlogram.  The absolute value of the correlogram can be useful for
%     finding delay times between records/signals that may or may not have 
%     opposite polarity.  A good case for this is when correlating a global
%     set of body wave phases that were recorded in different portions of 
%     the earthquake's focal sphere.  The unaltered case is useful for
%     times when you know there is no polarity difference, for example when
%     correlating teleseismic phases in a local array.

% check number of arguments
error(nargchk(1,9,nargin));

% defaults
if(nargin<9 || isempty(pow2pad)); pow2pad=1; end
if(nargin<8 || isempty(vectorxc)); vectorxc=false; end
if(nargin<7 || isempty(absxc)); absxc=false; end
if(nargin<6 || isempty(normxc)); normxc=false; end
if(nargin<5                   ); lags=[]; end
if(nargin<4 || isempty(spacing)); spacing=1; end
if(nargin<3 || isempty(npeaks)); npeaks=1; end

% checks
if(~isnumeric(pow2pad) || ~isscalar(pow2pad) || (fix(pow2pad)-pow2pad)~=0)
    error('SAClab:mcxc:badInput','POW2PAD must be an integer scalar'); 
end
if(~isnumeric(spacing) || ~isscalar(spacing) ...
        || spacing<1 || (fix(spacing)-spacing)~=0)
    error('SAClab:mcxc:badInput','SPACING must be an integer scalar >= 1'); 
end
if(~isnumeric(npeaks) || ~isscalar(npeaks) ...
        || npeaks<0 || (fix(npeaks)-npeaks)~=0)
    error('SAClab:mcxc:badInput','NPEAKS must be an integer scalar >= 0'); 
end
if(~isnumeric(lags) || (~isempty(lags) && ~isvector(lags))...
        || length(lags)>2 || any((fix(lags)-lags)~=0))
    error('SAClab:mcxc:badInput',...
        'LAGS must be an vector of integers with 0 to 2 elements');
end
if(~islogical(normxc))
    error('SAClab:mcxc:badInput','NORMXC must be logical');
end
if(~islogical(absxc))
    error('SAClab:mcxc:badInput','ABSXC must be logical'); 
end
if(~islogical(vectorxc))
    error('SAClab:mcxc:badInput','VECTORXC must be logical'); 
end
spacing=spacing-1;  % Intervals to samples
lags=sort(lags);    % Force order

% x or x+y
if(nargin==1 || isempty(y))
    [cg,lg,pg]=smcxc(x,npeaks,spacing,lags,normxc,absxc,vectorxc,pow2pad);
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
    Fx=fft(x,nFFT); Fy=fft(y,nFFT); CFy=conj(Fy);
    
    % get peaks
    if(npeaks)
        % vectorize
        if(vectorxc)
            % normalize
            if(normxc)
                % get autocorrelations
                ACFx=ifft(Fx.*conj(Fx));
                ACFy=ifft(Fy.*CFy);
                ZLACx=sqrt(ACFx(1,:));
                ZLACy=sqrt(ACFy(1,:));
                clear ACFx ACFy Fy
                
                % absolute peaks
                if(absxc)
                    [cg,lg,pg]=mcxcmncapv(Fx,CFy,ZLACx,ZLACy,nFFT,sx,sy,npeaks,spacing,lags);
                % positive peaks
                else
                    [cg,lg]=mcxcmncpv(Fx,CFy,ZLACx,ZLACy,nFFT,sx,sy,npeaks,spacing,lags);
                    if(nargout>2); pg=ones(sx(2),sy(2),npeaks); end
                end
            % dont normalize
            else
                clear Fy
                
                % absolute peaks
                if(absxc)
                    [cg,lg,pg]=mcxcmcapv(Fx,CFy,nFFT,sx,sy,npeaks,spacing,lags);
                % positive peaks
                else
                    [cg,lg]=mcxcmcpv(Fx,CFy,nFFT,sx,sy,npeaks,spacing,lags);
                    if(nargout>2); pg=ones(sx(2),sy(2),npeaks); end
                end
            end
        % loop
        else
            % normalize
            if(normxc)
                % get autocorrelations
                ACFx=ifft(Fx.*conj(Fx));
                ACFy=ifft(Fy.*CFy);
                ZLACx=sqrt(ACFx(1,:));
                ZLACy=sqrt(ACFy(1,:));
                clear ACFx ACFy Fy
                
                % absolute peaks
                if(absxc)
                    [cg,lg,pg]=mcxcmncap(Fx,CFy,ZLACx,ZLACy,nFFT,sx,sy,npeaks,spacing,lags);
                % positive peaks
                else
                    [cg,lg]=mcxcmncp(Fx,CFy,ZLACx,ZLACy,nFFT,sx,sy,npeaks,spacing,lags);
                    if(nargout>2); pg=ones(sx(2),sy(2),npeaks); end
                end
            % dont normalize
            else
                clear Fy
                
                % absolute peaks
                if(absxc)
                    [cg,lg,pg]=mcxcmcap(Fx,CFy,nFFT,sx,sy,npeaks,spacing,lags);
                % positive peaks
                else
                    [cg,lg]=mcxcmcp(Fx,CFy,nFFT,sx,sy,npeaks,spacing,lags);
                    if(nargout>2); pg=ones(sx(2),sy(2),npeaks); end
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
                % get autocorrelations
                ACFx=ifft(Fx.*conj(Fx));
                ACFy=ifft(Fy.*CFy);
                ZLACx=sqrt(ACFx(1,:));
                ZLACy=sqrt(ACFy(1,:));
                clear ACFx ACFy Fy
                [cg,lg]=mcxcncv(Fx,CFy,ZLACx,ZLACy,nFFT,sx,sy,lags);
            % dont normalize
            else
                clear Fy
                [cg,lg]=mcxccv(Fx,CFy,nFFT,sx,sy,lags);
            end
        % loop
        else
            % normalize
            if(normxc)
                % get autocorrelations
                ACFx=ifft(Fx.*conj(Fx));
                ACFy=ifft(Fy.*CFy);
                ZLACx=sqrt(ACFx(1,:));
                ZLACy=sqrt(ACFy(1,:));
                clear ACFx ACFy Fy
                [cg,lg]=mcxcnc(Fx,CFy,ZLACx,ZLACy,nFFT,sx,sy,lags);
            % dont normalize
            else
                clear Fy
                [cg,lg]=mcxcc(Fx,CFy,nFFT,sx,sy,lags);
            end
        end
    end
end

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
len=length(sh);

% preallocate output
cg=zeros(sx(2),sy(2),len);

end

function [cg,lags]=mcxcnc(F,CF,ZLACx,ZLACy,nFFT,sx,sy,lagrng)
%MCXCNC    mcxc 4 normalized correlograms - looped

% preallocate just about everything
[cg,lags,shift]=mcxcc_alloc(nFFT,sx,sy,lagrng);

% get correlograms
for i=1:sx(2)
    for j=1:sy(2)
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),'symmetric');
        cg(i,j,:)=XCF(shift)./(ZLACx(i)*ZLACy(j));
    end
end

end

function [cg,lags]=mcxcc(F,CF,nFFT,sx,sy,lagrng)
%MCXCC    mcxc 4 correlograms - looped

% preallocate just about everything
[cg,lags,shift]=mcxcc_alloc(nFFT,sx,sy,lagrng);

% get correlograms
for i=1:sx(2)
    for j=1:sy(2)
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),'symmetric');
        cg(i,j,:)=XCF(shift);
    end
end

end

function [cg,lags]=mcxcncv(F,CF,ZLACx,ZLACy,nFFT,sx,sy,lagrng)
%MCXCNCV    mcxc 4 normalized correlograms - vectorized

% allocate
[cg,lags,shift,len]=mcxcc_alloc(nFFT,sx,sy,lagrng);
ZLACy=ZLACy(ones(len,1),:);
cg=permute(cg,[1 3 2]);

% get correlograms
for i=1:sx(2)
    % do all possible pairings against record i
    XCF=ifft(F(:,i*ones(1,sy(2))).*CF,'symmetric');
    
    % put zero lag at center and normalize
    cg(i,:,:)=XCF(shift,:)./(ZLACx(i)*ZLACy);
end
cg=permute(cg,[1 3 2]);

end

function [cg,lags]=mcxccv(F,CF,nFFT,sx,sy,lagrng)
%MCXCCV    mcxc 4 normalized correlograms - vectorized

% allocate
[cg,lags,shift]=mcxcc_alloc(nFFT,sx,sy,lagrng);
cg=permute(cg,[1 3 2]);

% get correlograms
for i=1:sx(2)
    % do all possible pairings against record i
    XCF=ifft(F(:,i*ones(1,sy(2))).*CF,'symmetric');
    
    % put zero lag at center
    cg(i,:,:)=XCF(shift,:);
end
cg=permute(cg,[1 3 2]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                      FOR LOOP MCXC SECTION                      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cg,lg,lags,sh,len,range]=mcxc_alloc(nFFT,sx,sy,n,spacing,lagrng)
% allocation for mcxc - looped

% shifts and lags
sh=[(sx(1):-1:1)  (nFFT:-1:nFFT-sy(1)+2)];
lags=(-sx(1)+1:sy(1)-1).';
glags=(lags>=lagrng(1) & lags<=lagrng(2));
lags=lags(glags);
sh=sh(glags);
len=length(sh);
range=(-spacing:spacing).';

% preallocate output
cg=zeros(sx(2),sy(2),n); lg=cg;

end

function [cg,lg]=mcxcmncp(F,CF,ZLACx,ZLACy,nFFT,sx,sy,n,spacing,lagrng)
%MCXCMNCP    mcxc 4 multiple normalized correlogram peaks - looped

% preallocate just about everything
[cg,lg,lags,shift,len,range]=mcxc_alloc(nFFT,sx,sy,n,spacing,lagrng);

% get correlation peaks
for i=1:sx(2)
    for j=1:sy(2)
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),'symmetric');
        XCFT=XCF(shift)./(ZLACx(i)*ZLACy(j));
        
        % pick peak
        [cg(i,j,1),imxc]=max(XCFT);
        lg(i,j,1)=lags(imxc);
        
        % multi-peak picker
        for k=2:n
            % zero previous peak
            izero=imxc+range;
            XCFT(izero(izero>0 & izero<=len))=0;
            
            % pick peak
            [cg(i,j,k),imxc]=max(XCFT);
            lg(i,j,k)=lags(imxc);
        end
    end
end

end

function [cg,lg]=mcxcmcp(F,CF,nFFT,sx,sy,n,spacing,lagrng)
%MCXCMCP    mcxc 4 multiple correlogram peaks - looped

% preallocate just about everything
[cg,lg,lags,shift,len,range]=mcxc_alloc(nFFT,sx,sy,n,spacing,lagrng);

% get correlation peaks
for i=1:sx(2)
    for j=1:sy(2)
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),'symmetric');
        XCFT=XCF(shift);
        
        % pick peak
        [cg(i,j,1),imxc]=max(XCFT);
        lg(i,j,1)=lags(imxc);
        
        % multi-peak picker
        for k=2:n
            % zero previous peak
            izero=imxc+range;
            XCFT(izero(izero>0 & izero<=len))=0;
            
            % pick peak
            [cg(i,j,k),imxc]=max(XCFT);
            lg(i,j,k)=lags(imxc);
        end
    end
end

end

function [cg,lg,pg]=mcxcmncap(F,CF,ZLACx,ZLACy,nFFT,sx,sy,n,spacing,lagrng)
%MCXCMNCAP    mcxc 4 multiple normalized correlogram absolute peaks - looped

% preallocate just about everything
[cg,lg,lags,shift,len,range]=mcxc_alloc(nFFT,sx,sy,n,spacing,lagrng);
pg=ones(sx(2),sy(2),n);

% get correlation peaks
for i=1:sx(2)
    for j=1:sy(2)
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),'symmetric');
        XCFT=XCF(shift)./(ZLACx(i)*ZLACy(j));
        AXCFT=abs(XCFT);
        
        % pick peak
        [cg(i,j,1),imxc]=max(AXCFT);
        lg(i,j,1)=lags(imxc);
        pg(i,j,1)=cg(i,j,1)/XCFT(imxc);
        
        % multi-peak picker
        for k=2:n
            % zero previous peak
            izero=imxc+range;
            AXCFT(izero(izero>0 & izero<=len))=0;
            
            % pick peak
            [cg(i,j,k),imxc]=max(AXCFT);
            lg(i,j,k)=lags(imxc);
            pg(i,j,k)=cg(i,j,k)/XCFT(imxc);
        end
    end
end

end

function [cg,lg,pg]=mcxcmcap(F,CF,nFFT,sx,sy,n,spacing,lagrng)
%MCXCMCAP    mcxc 4 multiple correlogram absolute peaks - looped

% preallocate just about everything
[cg,lg,lags,shift,len,range]=mcxc_alloc(nFFT,sx,sy,n,spacing,lagrng);
pg=ones(sx(2),sy(2),n);

% get correlation peaks
for i=1:sx(2)
    for j=1:sy(2)
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),'symmetric');
        XCFT=XCF(shift);
        AXCFT=abs(XCFT);
        
        % pick peak
        [cg(i,j,1),imxc]=max(AXCFT);
        lg(i,j,1)=lags(imxc);
        pg(i,j,1)=cg(i,j,1)/XCFT(imxc);
        
        % multi-peak picker
        for k=2:n
            % zero previous peak
            izero=imxc+range;
            AXCFT(izero(izero>0 & izero<=len))=0;
            
            % pick peak
            [cg(i,j,k),imxc]=max(AXCFT);
            lg(i,j,k)=lags(imxc);
            pg(i,j,k)=cg(i,j,k)/XCFT(imxc);
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                      VECTORIZED MCXC SECTION                    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cg,lg,lags,sh,len,range,per21,cols]=mcxcv_alloc(nFFT,sx,sy,n,spacing,lagrng)
% allocation for mcxc - vectorized

% shifts and lags
sh=[(sx(1):-1:1)  (nFFT:-1:nFFT-sy(1)+2)];
lags=(-sx(1)+1:sy(1)-1).';
glags=(lags>=lagrng(1) & lags<=lagrng(2));
lags=lags(glags);
sh=sh(glags);
len=length(sh);
range=(-spacing:spacing).';
per21=2*spacing+1;

% preallocate output
cg=zeros(sx(2),sy(2),n); lg=cg;

% preallocate zeroer
cols=repmat(1:sy(2),per21,1);

end

function [cg,lg]=mcxcmncpv(F,CF,ZLACx,ZLACy,nFFT,sx,sy,n,spacing,lagrng)
%MCXCMNCPV    mcxc 4 multiple normalized correlogram peaks - vectorized

% allocate
[cg,lg,lags,shift,len,range,per21,cols]=mcxcv_alloc(nFFT,sx,sy,n,spacing,lagrng);
ZLACy=ZLACy(ones(len,1),:);

% get correlation peaks
for i=1:sx(2)
    % do all possible pairings against record i
    XCF=ifft(F(:,i*ones(1,sy(2))).*CF,'symmetric');
    
    % put zero lag at center and normalize
    XCFT=XCF(shift,:)./(ZLACx(i)*ZLACy);
    
    % find peak
    [cg(i,:,1),imxc]=max(XCFT);
    lg(i,:,1)=lags(imxc);
    
    % multi-peak picker
    for j=2:n
        % zero out peak
        rows=imxc(ones(per21,1),:)+range(:,ones(1,sy(2)));
        good=(rows>0 & rows<=len);
        XCFT(sub2ind([len sy(2)],rows(good),cols(good)))=0;
        
        % find peak
        [cg(i,:,j),imxc]=max(XCFT);
        lg(i,:,j)=lags(imxc);
    end
end

end

function [cg,lg]=mcxcmcpv(F,CF,nFFT,sx,sy,n,spacing,lagrng)
%MCXCMCPV    mcxc 4 multiple correlogram peaks - vectorized

% allocate
[cg,lg,lags,shift,len,range,per21,cols]=mcxcv_alloc(nFFT,sx,sy,n,spacing,lagrng);

% get correlation peaks
for i=1:sx(2)
    % do all possible pairings against record i
    XCF=ifft(F(:,i*ones(1,sy(2))).*CF,'symmetric');
    
    % put zero lag at center
    XCFT=XCF(shift,:);
    
    % find peak
    [cg(i,:,1),imxc]=max(XCFT);
    lg(i,:,1)=lags(imxc);
    
    % multi-peak picker
    for j=2:n
        % zero out peak
        rows=imxc(ones(per21,1),:)+range(:,ones(1,sy(2)));
        good=(rows>0 & rows<=len);
        XCFT(sub2ind([len sy(2)],rows(good),cols(good)))=0;
        
        % find peak
        [cg(i,:,j),imxc]=max(XCFT);
        lg(i,:,j)=lags(imxc);
    end
end

end

function [cg,lg,pg]=mcxcmncapv(F,CF,ZLACx,ZLACy,nFFT,sx,sy,n,spacing,lagrng)
%MCXCMNCAPV    mcxc 4 multiple normalized correlogram absolute peaks - vectorized

% allocate
[cg,lg,lags,shift,len,range,per21,cols]=mcxcv_alloc(nFFT,sx,sy,n,spacing,lagrng);
ZLACy=ZLACy(ones(len,1),:);
pg=ones(sx(2),sy(2),n);

% get correlation peaks
for i=1:sx(2)
    % do all possible pairings against record i
    XCF=ifft(F(:,i*ones(1,sy(2))).*CF,'symmetric');
    
    % put zero lag at center and make real
    XCFT=XCF(shift,:)./(ZLACx(i)*ZLACy);
    
    % also find absolute values
    AXCFT=abs(XCFT);
    
    % find peak
    [cg(i,:,1),imxc]=max(AXCFT);
    lg(i,:,1)=lags(imxc);
    pg(i,:,1)=cg(i,:,1)./XCFT((0:sy(2)-1)*len+imxc);
    
    % multi-peak picker
    for j=2:n
        % zero out peak
        rows=imxc(ones(per21,1),:)+range(:,ones(1,sy(2)));
        good=(rows>0 & rows<=len);
        AXCFT(sub2ind([len sy(2)],rows(good),cols(good)))=0;
        
        % find peak
        [cg(i,:,j),imxc]=max(AXCFT);
        lg(i,:,j)=lags(imxc);
        pg(i,:,j)=cg(i,:,j)./XCFT((0:sy(2)-1)*len+imxc);
    end
end

end

function [cg,lg,pg]=mcxcmcapv(F,CF,nFFT,sx,sy,n,spacing,lagrng)
%MCXCMCAPV    mcxc 4 multiple correlogram absolute peaks - vectorized

% allocate
[cg,lg,lags,shift,len,range,per21,cols]=mcxcv_alloc(nFFT,sx,sy,n,spacing,lagrng);
pg=ones(sx(2),sy(2),n);

% get correlation peaks
for i=1:sx(2)
    % do all possible pairings against record i
    XCF=ifft(F(:,i*ones(1,sy(2))).*CF,'symmetric');
    
    % put zero lag at center and make real
    XCFT=XCF(shift,:);
    
    % also find absolute values
    AXCFT=abs(XCFT);
    
    % find peak
    [cg(i,:,1),imxc]=max(AXCFT);
    lg(i,:,1)=lags(imxc);
    pg(i,:,1)=cg(i,:,1)./XCFT((0:sy(2)-1)*len+imxc);
    
    % multi-peak picker
    for j=2:n
        % zero out peak
        rows=imxc(ones(per21,1),:)+range(:,ones(1,sy(2)));
        good=(rows>0 & rows<=len);
        AXCFT(sub2ind([len sy(2)],rows(good),cols(good)))=0;
        
        % find peak
        [cg(i,:,j),imxc]=max(AXCFT);
        lg(i,:,j)=lags(imxc);
        pg(i,:,j)=cg(i,:,j)./XCFT((0:sy(2)-1)*len+imxc);
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

function [cv,lv,pv]=smcxc(x,npeaks,spacing,lags,normxc,absxc,vectorxc,pow2pad)
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
CF=conj(F);

% get peaks
if(npeaks)
    % vectorize
    if(vectorxc)
        % normalize
        if(normxc)
            % get autocorrelations
            ACF=ifft(F.*CF,'symmetric');
            ZLAC=sqrt(ACF(1,:)); clear ACF
            
            % absolute peaks
            if(absxc)
                [cv,lv,pv]=smcxcmncapv(F,CF,ZLAC,nFFT,sx,npeaks,spacing,lags);
            % positive peaks
            else
                [cv,lv]=smcxcmncpv(F,CF,ZLAC,nFFT,sx,npeaks,spacing,lags);
                if(nargout>2); pv=ones(txc,1,npeaks); end
            end
        % dont normalize
        else
            % absolute peaks
            if(absxc)
                [cv,lv,pv]=smcxcmcapv(F,CF,nFFT,sx,npeaks,spacing,lags);
            % positive peaks
            else
                [cv,lv]=smcxcmcpv(F,CF,nFFT,sx,npeaks,spacing,lags);
                if(nargout>2); pv=ones(txc,1,npeaks); end
            end
        end
    % loop
    else
        % normalize
        if(normxc)
            % get autocorrelations
            ACF=ifft(F.*CF,'symmetric');
            ZLAC=sqrt(ACF(1,:)); clear ACF
            
            % absolute peaks
            if(absxc)
                [cv,lv,pv]=smcxcmncap(F,CF,ZLAC,nFFT,sx,npeaks,spacing,lags);
            % positive peaks
            else
                [cv,lv]=smcxcmncp(F,CF,ZLAC,nFFT,sx,npeaks,spacing,lags);
                if(nargout>2); pv=ones(txc,1,npeaks); end
            end
        % dont normalize
        else
            % absolute peaks
            if(absxc)
                [cv,lv,pv]=smcxcmcap(F,CF,nFFT,sx,npeaks,spacing,lags);
            % positive peaks
            else
                [cv,lv]=smcxcmcp(F,CF,nFFT,sx,npeaks,spacing,lags);
                if(nargout>2); pv=ones(txc,1,npeaks); end
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
            % get autocorrelations
            ACF=ifft(F.*CF,'symmetric');
            ZLAC=sqrt(ACF(1,:)); clear ACF
            [cv,lv]=smcxcncv(F,CF,ZLAC,nFFT,sx,lags);
        % dont normalize
        else
            [cv,lv]=smcxccv(F,CF,nFFT,sx,lags);
        end
    % loop
    else
        % normalize
        if(normxc)
            % get autocorrelations
            ACF=ifft(F.*CF,'symmetric');
            ZLAC=sqrt(ACF(1,:)); clear ACF
            [cv,lv]=smcxcnc(F,CF,ZLAC,nFFT,sx,lags);
        % dont normalize
        else
            [cv,lv]=smcxcc(F,CF,nFFT,sx,lags);
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 CORRELOGRAM ONLY SMCXC SECTION                  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cv,lags,sh,c]=smcxcc_alloc(nFFT,sx,lagrng)
% allocation for smcxcc - looped

% shifts and lags
sh=[(sx(1):-1:1)  (nFFT:-1:nFFT-sx(1)+2)];
lags=(-sx(1)+1:sx(1)-1).';
glags=(lags>=lagrng(1) & lags<=lagrng(2));
lags=lags(glags);
sh=sh(glags);
len=length(sh);
ilen=sx(2)-1:-1:1;
c=[0 cumsum(ilen)];
txc=(sx(2)^2-sx(2))/2;

% preallocate output
cv=zeros(txc,1,len);

end

function [cv,lags]=smcxcnc(F,CF,ZLAC,nFFT,sx,lagrng)
%SMCXCNC    smcxc 4 normalized correlograms - looped

% preallocate just about everything
[cv,lags,shift,c]=smcxcc_alloc(nFFT,sx,lagrng);

% get correlograms
for i=1:sx(2)-1
    for j=i+1:sx(2)
        k=c(i)+j-i;
        
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),'symmetric');
        cv(k,1,:)=XCF(shift)./(ZLAC(i)*ZLAC(j));
    end
end

end

function [cv,lags]=smcxcc(F,CF,nFFT,sx,lagrng)
%SMCXCC    smcxc 4 correlograms - looped

% preallocate just about everything
[cv,lags,shift,c]=smcxcc_alloc(nFFT,sx,lagrng);

% get correlograms
for i=1:sx(2)-1
    for j=i+1:sx(2)
        k=c(i)+j-i;
        
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),'symmetric');
        cv(k,1,:)=XCF(shift);
    end
end

end

function [cv,XCF,lags,sh,len,ilen,c,c2]=smcxccv_alloc(nFFT,sx,lagrng)
% allocation for smcxcc - vectorized

% preallocate shifts and lags
sh=[(sx(1):-1:1)  (nFFT:-1:nFFT-sx(1)+2)];
lags=(-sx(1)+1:sx(1)-1).';
glags=(lags>=lagrng(1) & lags<=lagrng(2));
lags=lags(glags);
sh=sh(glags);
len=length(sh);
ilen=sx(2)-1:-1:1;
c2=cumsum(ilen);
c=[1 c2+1];
txc=(sx(2)^2-sx(2))/2;

% preallocate output
cv=zeros(txc,1,len);

% preallocate correlograms
XCF=zeros(nFFT,sx(2));

end

function [cv,lags]=smcxcncv(F,CF,ZLAC,nFFT,sx,lagrng)
%SMCXCNCV    smcxc 4 normalized correlograms - vectorized

% preallocate just about everything
[cv,XCF,lags,shift,len,ilen,c,c2]=smcxccv_alloc(nFFT,sx,lagrng);
ZLAC2=ZLAC(ones(len,1),:);
cv=permute(cv,[3 2 1]);

% get correlograms
for i=1:sx(2)-1
    % indicies
    i2=i+1:sx(2);
    
    % do all possible pairings against record i
    XCF(:,i2)=ifft(F(:,i*ones(1,ilen(i))).*CF(:,i2),'symmetric');
    
    % put zero lag at center & normalize
    cv(:,1,c(i):c2(i))=XCF(shift,i2)./(ZLAC(i)*ZLAC2(:,i2));
end
cv=permute(cv,[3 2 1]);

end

function [cv,lags]=smcxccv(F,CF,nFFT,sx,lagrng)
%SMCXCCV    smcxc 4 correlograms - vectorized

% preallocate just about everything
[cv,XCF,lags,shift,len,ilen,c,c2]=smcxccv_alloc(nFFT,sx,lagrng);
cv=permute(cv,[3 2 1]);

% get correlograms
for i=1:sx(2)-1
    % indicies
    i2=i+1:sx(2);
    
    % do all possible pairings against record i
    XCF(:,i2)=ifft(F(:,i*ones(1,ilen(i))).*CF(:,i2),'symmetric');
    
    % put zero lag at center
    cv(:,1,c(i):c2(i))=XCF(shift,i2);
end
cv=permute(cv,[3 2 1]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                       LOOPED SMCXC SECTION                      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cv,lv,lags,sh,len,range,per21,c,txc]=smcxc_alloc(nFFT,sx,n,spacing,lagrng)
% allocation for smcxc - looped

% preallocate shifts and lags
sh=[(sx(1):-1:1)  (nFFT:-1:nFFT-sx(1)+2)];
lags=(-sx(1)+1:sx(1)-1).';
glags=(lags>=lagrng(1) & lags<=lagrng(2));
lags=lags(glags);
sh=sh(glags);
len=length(sh);
range=(-spacing:spacing).';
per21=2*spacing+1;
ilen=sx(2)-1:-1:1;
c=[0 cumsum(ilen)];
txc=(sx(2)^2-sx(2))/2;

% preallocate output
cv=zeros(txc,1,n); lv=cv;

end

function [cv,lv]=smcxcmncp(F,CF,ZLAC,nFFT,sx,n,spacing,lagrng)
%SMCXCMNCP    smcxc 4 multiple normalized correlogram peaks - looped

% preallocate just about everything
[cv,lv,lags,shift,len,range,per21,c]=smcxc_alloc(nFFT,sx,n,spacing,lagrng);

% get correlation peaks
for i=1:sx(2)-1
    for j=i+1:sx(2)
        k=c(i)+j-i;
        
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),'symmetric');
        XCFT=XCF(shift)./(ZLAC(i)*ZLAC(j));
        
        % pick peak
        [cv(k,1,1),imxc]=max(XCFT);
        lv(k,1,1)=lags(imxc);
        
        % multi-peak picker
        for p=2:n
            % zero previous peak
            izero=imxc+range;
            XCFT(izero(izero>0 & izero<=len))=0;
            
            % pick peak
            [cv(k,1,p),imxc]=max(XCFT);
            lv(k,1,p)=lags(imxc);
        end
    end
end

end

function [cv,lv]=smcxcmcp(F,CF,nFFT,sx,n,spacing,lagrng)
%SMCXCMCP    smcxc 4 multiple correlogram peaks - looped

% preallocate just about everything
[cv,lv,lags,shift,len,range,per21,c]=smcxc_alloc(nFFT,sx,n,spacing,lagrng);

% get correlation peaks
for i=1:sx(2)-1
    for j=i+1:sx(2)
        k=c(i)+j-i;
        
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),'symmetric');
        XCFT=XCF(shift);
        
        % pick peak
        [cv(k,1,1),imxc]=max(XCFT);
        lv(k,1,1)=lags(imxc);
        
        % multi-peak picker
        for p=2:n
            % zero previous peak
            izero=imxc+range;
            XCFT(izero(izero>0 & izero<=len))=0;
            
            % pick peak
            [cv(k,1,p),imxc]=max(XCFT);
            lv(k,1,p)=lags(imxc);
        end
    end
end

end

function [cv,lv,pv]=smcxcmncap(F,CF,ZLAC,nFFT,sx,n,spacing,lagrng)
%SMCXCMNCAP    smcxc 4 multiple normalized correlogram absolute peaks - looped

% preallocate just about everything
[cv,lv,lags,shift,len,range,per21,c,txc]=smcxc_alloc(nFFT,sx,n,spacing,lagrng);
pv=ones(txc,1,n);

% get correlation peaks
for i=1:sx(2)-1
    for j=i+1:sx(2)
        k=c(i)+j-i;
        
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),'symmetric');
        XCFT=XCF(shift)./(ZLAC(i)*ZLAC(j));
        AXCFT=abs(XCFT);
        
        % pick peak
        [cv(k,1,1),imxc]=max(AXCFT);
        lv(k,1,1)=lags(imxc);
        pv(k,1,1)=cv(k,1,1)/XCFT(imxc);
        
        % multi-peak picker
        for p=2:n
            % zero previous peak
            izero=imxc+range;
            AXCFT(izero(izero>0 & izero<=len))=0;
            
            % pick peak
            [cv(k,1,p),imxc]=max(AXCFT);
            lv(k,1,p)=lags(imxc);
            pv(k,1,p)=cv(k,1,p)/XCFT(imxc);
        end
    end
end

end

function [cv,lv,pv]=smcxcmcap(F,CF,nFFT,sx,n,spacing,lagrng)
%SMCXCMCAP    smcxc 4 multiple correlogram absolute peaks - looped

% preallocate just about everything
[cv,lv,lags,shift,len,range,per21,c,txc]=smcxc_alloc(nFFT,sx,n,spacing,lagrng);
pv=ones(txc,1,n);

% get correlation peaks
for i=1:sx(2)-1
    for j=i+1:sx(2)
        k=c(i)+j-i;
        
        % correlate
        XCF=ifft(F(:,i).*CF(:,j),'symmetric');
        XCFT=XCF(shift);
        AXCFT=abs(XCFT);
        
        % pick peak
        [cv(k,1,1),imxc]=max(AXCFT);
        lv(k,1,1)=lags(imxc);
        pv(k,1,1)=cv(k,1,1)/XCFT(imxc);
        
        % multi-peak picker
        for p=2:n
            % zero previous peak
            izero=imxc+range;
            AXCFT(izero(izero>0 & izero<=len))=0;
            
            % pick peak
            [cv(k,1,p),imxc]=max(AXCFT);
            lv(k,1,p)=lags(imxc);
            pv(k,1,p)=cv(k,1,p)/XCFT(imxc);
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                    VECTORIZED SMCXC SECTION                     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cv,lv,XCF,XCFT,lags,sh,len,range,per21,ilen,c,c2,imxc,rows,cols,txc]=smcxcv_alloc(nFFT,sx,n,spacing,lagrng)
% allocation for smcxc - vectorized

% preallocate shifts and lags
sh=[(sx(1):-1:1)  (nFFT:-1:nFFT-sx(1)+2)];
lags=(-sx(1)+1:sx(1)-1).';
glags=(lags>=lagrng(1) & lags<=lagrng(2));
lags=lags(glags);
sh=sh(glags);
len=length(sh);
range=(-spacing:spacing).';
per21=2*spacing+1;
ilen=sx(2)-1:-1:1;
c2=cumsum(ilen);
c=[1 c2+1];
txc=(sx(2)^2-sx(2))/2;
imxc=zeros(1,sx(2));

% preallocate output
cv=zeros(txc,1,n); lv=cv;

% preallocate correlograms
XCF=zeros(nFFT,sx(2));
XCFT=zeros(len,sx(2));

% preallocate peak zeroer
rows=zeros(per21,sx(2));
cols=repmat(1:sx(2),per21,1);

end

function [cv,lv]=smcxcmncpv(F,CF,ZLAC,nFFT,sx,n,spacing,lagrng)
%SMCXCMNCPV    smcxc 4 multiple normalized correlogram peaks - vectorized

% preallocate just about everything
[cv,lv,XCF,XCFT,lags,shift,len,range,per21,ilen,c,c2,imxc,rows,cols]=...
    smcxcv_alloc(nFFT,sx,n,spacing,lagrng);
ZLAC2=ZLAC(ones(len,1),:);

% get correlation peaks
for i=1:sx(2)-1
    % indicies
    i2=i+1:sx(2);
    
    % do all possible pairings against record i
    XCF(:,i2)=ifft(F(:,i*ones(1,ilen(i))).*CF(:,i2),'symmetric');
    
    % put zero lag at center & normalize
    XCFT(:,i2)=XCF(shift,i2)./(ZLAC(i)*ZLAC2(:,i2));
    
    % find peak (for each pairing)
    [cv(c(i):c2(i),1,1),imxc(i2)]=max(XCFT(:,i2));
    lv(c(i):c2(i),1,1)=lags(imxc(i2));
    
    % multi-peak picker
    for j=2:n
        % zero out peak
        rows(:,i2)=imxc(ones(per21,1),i2)+range(:,ones(1,ilen(i)));
        good=(rows>0 & rows<=len);
        XCFT(sub2ind([len sx(2)],rows(good),cols(good)))=0;
        
        % find next peak
        [cv(c(i):c2(i),1,j),imxc(i2)]=max(XCFT(:,i2));
        lv(c(i):c2(i),1,j)=lags(imxc(i2));
    end
end

end

function [cv,lv]=smcxcmcpv(F,CF,nFFT,sx,n,spacing,lagrng)
%SMCXCMCPV    smcxc 4 multiple correlogram peaks - vectorized

% preallocate just about everything
[cv,lv,XCF,XCFT,lags,shift,len,range,per21,ilen,c,c2,imxc,rows,cols]=...
    smcxcv_alloc(nFFT,sx,n,spacing,lagrng);

% get correlation peaks
for i=1:sx(2)-1
    % indicies
    i2=i+1:sx(2);
    
    % do all possible pairings against record i
    XCF(:,i2)=ifft(F(:,i*ones(1,ilen(i))).*CF(:,i2),'symmetric');
    
    % put zero lag at center
    XCFT(:,i2)=XCF(shift,i2);
    
    % find peak (for each pairing)
    [cv(c(i):c2(i),1,1),imxc(i2)]=max(XCFT(:,i2));
    lv(c(i):c2(i),1,1)=lags(imxc(i2));
    
    % multi-peak picker
    for j=2:n
        % zero out peak
        rows(:,i2)=imxc(ones(per21,1),i2)+range(:,ones(1,ilen(i)));
        good=(rows>0 & rows<=len);
        XCFT(sub2ind([len sx(2)],rows(good),cols(good)))=0;
        
        % find next peak
        [cv(c(i):c2(i),1,j),imxc(i2)]=max(XCFT(:,i2));
        lv(c(i):c2(i),1,j)=lags(imxc(i2));
    end
end

end

function [cv,lv,pv]=smcxcmncapv(F,CF,ZLAC,nFFT,sx,n,spacing,lagrng)
%SMCXCMNCAPV    smcxc 4 multiple normalized correlogram absolute peaks - vectorized

% preallocate just about everything
[cv,lv,XCF,XCFT,lags,shift,len,range,per21,ilen,c,c2,imxc,rows,cols,txc]=...
    smcxcv_alloc(nFFT,sx,n,spacing,lagrng);
ZLAC2=ZLAC(ones(len,1),:);
pv=ones(txc,1,n);
AXCFT=XCFT;

% get correlation peaks
for i=1:sx(2)-1
    % indicies
    i2=i+1:sx(2);
    
    % do all possible pairings against record i
    XCF(:,i2)=ifft(F(:,i*ones(1,ilen(i))).*CF(:,i2),'symmetric');
    
    % put zero lag at center & normalize
    XCFT(:,i2)=XCF(shift,i2)./(ZLAC(i)*ZLAC2(:,i2));
    
    % also find absolute values
    AXCFT(:,i2)=abs(XCFT(:,i2));
    
    % find peak (for each pairing)
    [cv(c(i):c2(i),1,1),imxc(i2)]=max(AXCFT(:,i2));
    lv(c(i):c2(i),1,1)=lags(imxc(i2));
    pv(c(i):c2(i),1,1)=cv(c(i):c2(i),1,1)./XCFT(imxc(i2)+(i2-1)*len).';
    
    % multi-peak picker
    for j=2:n
        % zero out peak
        rows(:,i2)=imxc(ones(per21,1),i2)+range(:,ones(1,ilen(i)));
        good=(rows>0 & rows<=len);
        AXCFT(sub2ind([len sx(2)],rows(good),cols(good)))=0;
        
        % find next peak
        [cv(c(i):c2(i),1,j),imxc(i2)]=max(AXCFT(:,i2));
        lv(c(i):c2(i),1,j)=lags(imxc(i2));
        pv(c(i):c2(i),1,j)=cv(c(i):c2(i),1,j)./XCFT(imxc(i2)+(i2-1)*len).';
    end
end

end

function [cv,lv,pv]=smcxcmcapv(F,CF,nFFT,sx,n,spacing,lagrng)
%SMCXCMCAPV    smcxc 4 multiple correlogram absolute peaks - vectorized

% preallocate just about everything
[cv,lv,XCF,XCFT,lags,shift,len,range,per21,ilen,c,c2,imxc,rows,cols,txc]=...
    smcxcv_alloc(nFFT,sx,n,spacing,lagrng);
pv=ones(txc,1,n);
AXCFT=XCFT;

% get correlation peaks
for i=1:sx(2)-1
    % indicies
    i2=i+1:sx(2);
    
    % do all possible pairings against record i
    XCF(:,i2)=ifft(F(:,i*ones(1,ilen(i))).*CF(:,i2),'symmetric');
    
    % put zero lag at center & normalize
    XCFT(:,i2)=XCF(shift,i2);
    
    % also find absolute values
    AXCFT(:,i2)=abs(XCFT(:,i2));
    
    % find peak (for each pairing)
    [cv(c(i):c2(i),1,1),imxc(i2)]=max(AXCFT(:,i2));
    lv(c(i):c2(i),1,1)=lags(imxc(i2));
    pv(c(i):c2(i),1,1)=cv(c(i):c2(i),1,1)./XCFT(imxc(i2)+(i2-1)*len).';
    
    % multi-peak picker
    for j=2:n
        % zero out peak
        rows(:,i2)=imxc(ones(per21,1),i2)+range(:,ones(1,ilen(i)));
        good=(rows>0 & rows<=len);
        AXCFT(sub2ind([len sx(2)],rows(good),cols(good)))=0;
        
        % find next peak
        [cv(c(i):c2(i),1,j),imxc(i2)]=max(AXCFT(:,i2));
        lv(c(i):c2(i),1,j)=lags(imxc(i2));
        pv(c(i):c2(i),1,j)=cv(c(i):c2(i),1,j)./XCFT(imxc(i2)+(i2-1)*len).';
    end
end

end
