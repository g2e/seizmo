function [data,fs,nyq]=iirfilter(data,type,style,corners,order,passes,ripple)
%IIRFILTER    Apply an IIR filter to SAClab data records
%
%    Description: Utilizing the set of supplied parameters IIR filters are
%     built to the desired specs and implemented on the data records.
%
%     Parameters are as follows:
%
%     Type:    'high', 'low', 'notch' or 'bandpass'
%
%     Style:   'butter', 'cheby1', 'cheby2' or 'ellip'
%
%     Corners: Array of pass corner(s) and optionally stop corner(s) - for
%              notch type filtering the opposite is true (and as such, all
%              of the rest of the explanation that follows).  If stop
%              corner(s) are given, then the order parameter is ignored 
%              and automatic filter order determination is done utilizing
%              those corners and the ripple parameters (if given).  If no
%              stop corner(s) are given and automatic order calculation is 
%              necessary (for instance if the ripple parameters are set), 
%              default stop corner(s) are set up.  These are set at 1/3 the 
%              passband width from the passband corner(s).  If this puts
%              the corners outside the bandwidth of the signal (0Hz to 
%              nyquist) the corner is adjusted to the nyquist/0Hz mark.  
%
%              For lowpass or highpass filtering 1 or 2 corners can be
%              defined.  Bandpass or notch filters require 2 or 4 corners.
%
%         ***  UNITS ARE IN Hz!!!  ***
%
%              Special Note: Due to the nature of its design, FOR THE 
%                            'CHEBY2' FILTER, THE STOPBAND CORNERS ARE 
%                            REQUIRED WHILE THE PASS-BAND CORNERS ARE 
%                            OPTIONAL.  Defaults when no passband corners 
%                            are given are set so that the transition width
%                            is 1/3 the passband width (to match the other
%                            filter style defaults).  Things get rather
%                            confusing for a notch cheby2 filter...
%
%                            Passband corners set the frequency at which
%                            the attenuation hits 3dB, while stopband
%                            corners set the 30dB position.  For a notch
%                            filter this is reversed.  The dB levels that
%                            correspond to the corners can be adjusted
%                            utilizing the ripple parameter.
%                          
%     Order:  Optional.  If undefined or 0, automatic filter order
%             calculations are made (utilizes the supplied/default corners
%             and ripple parameter to define a proper order).  Filter 
%             order is analogous to the sharpness/steepness/slope of the
%             transition from the passband to the stopband.  Higher orders 
%             provide better frequency resolution at the cost time 
%             resolution (aka ringing).
%
%     Passes: Optional.  Defaults to 1.  Accepts 1,2,3 or 4.  Options 3 is
%             for backwards filtering (option 4 forward filters after
%             backwards filtering - aka a 2-pass in reverse order).  Two
%             pass filtering is convenient for zero-phase filtering as
%             single pass IIR filtering causes phase dispersion.  Note that
%             2-pass filters make the signal acausal and should not be done
%             in conjunction with phase onset picking.
%
%     Ripple: Optional.  1 or 2 element array of ripple/attenuation values
%             indicating the maximum permissible passband loss and minimum 
%             stopband loss. Defaults are 3dB for the passband and 30dB for 
%             the stopband. Both values are used in automatic filter order 
%             determination.
%
%             Automatic filter order calculation will be done if both 
%             values are supplied unless the filter is elliptic. If any
%             ripple parameter is given for butter type, automatic filter
%             order calculation will be done.
%
%             For elliptic and cheby1 filtering if one value is given it is
%             taken as the passband ripple.  For cheby2 and butter it is
%             taken as the stopband ripple/attenuation.
%
%        ***  UNITS IN dB!!!  ***
%
%     fs:     Optional output. Filter structure(s).  Useful for visualizing 
%             the filter with the fvtool.  Remember that the filter is tied
%             to the sampling frequency of the records!
%
%     nyq:    Optional output. Nyquist frequencies for filter structure(s).
%
%
%    Usage: [data,fs]=iirfilter(data,type,style,corners,order,passes,ripple)
%
%
%    Examples:
%      4th order lowpass butter filter with a passband corner of 10s
%     [data]=sacfilter(data,'low','butter',1/10,4)
%
%      bandpass elliptic filter with corners of 15s and 12s using
%      default stop corners and ripples to determine filter order
%     [data]=sacfilter(data,'bandpass','ellip',[1/15 1/12])
%
%      2 pass bandpass elliptic filter with stop corners 
%      defined forcing order calc. 8th order is ignored but the
%      supplied 3dB passband ripple is included in the order calc.
%     [data]=sacfilter(data,'bandpass','ellip',[1/22 1/20 1/15 1/14],8,2,3)
%
%
%    Notes: For narrow band filters, PARTICULARLY AT LOW Hz, consider 
%           decimating the signal prior to filtering to aid in the accuracy
%           of the filter design.  This is because the analog to digital 
%           transformation of the filter design will lose accuracy if the 
%           features of the filter are tighter than the digitization level.
%           Basically try to keep the filter width above a few percent of 
%           the total bandwidth of the signal.  
%
%           Remember, length of the trace does not affect filter design,
%           but the sampling frequency does.  
%
%           Edge effect limiting in filtering is done by 'mirroring' the 
%           data record before filtering which serves as a pseudo-initial
%           condition determination.  The cost is that of filtering a
%           record twice as long.
%
%           This program will filter on any filetype.  The output for 
%           unevenly spaced or xyz files will not be accurate though!
%
%    See also: intrpol8, syncsr, deci, stretch, fourier, ifourier

% check number of inputs
error(nargchk(2,7,nargin))

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head') || ...
        ~isfield(data,'x'))
    error('data structure does not have proper fields')
end

% check filter
if(~any(strcmpi(style,{'butter' 'cheby1' 'cheby2' 'ellip'})) || ...
   ~any(strcmpi(type,{'low' 'high' 'notch' 'stop' 'bandpass'})))
        error('unknown filter')
end

% built-in defaults
auto=0;     % auto order det flag
tranbw=1/3; % trans/passband width ratio
mirror=1;   % mirror flag
pr=3;       % passband attenuation (dB)
sr=30;      % stopband attenuation (dB)
dp=1;       % pass option

% make delta groups
delta=gh(data,'delta');
gdt=unique(delta); ng=length(gdt); gi=cell(ng,1); nyq=1./(2*gdt);
for i=1:ng; gi{i}=find(delta==gdt(i)).'; end

% notch to stop
if(strcmpi(type,'notch')); type='stop'; end

% check order
if(nargin==4 || isempty(order) || order==0); auto=1; end

% check passes
if(nargin<6 || isempty(passes)); passes=dp;
elseif (~any(passes==[1 2 3 4])); error('1 or 2 passes only'); end

% check ripple
if(nargin==7 && ~isempty(ripple))
    ripple=sort(ripple(:).');
    nr=length(ripple);

    % check for negative
    if(any(ripple<=0))
        error('negative or zero ripple not allowed')
    end

    % assigning ripple parameters
    if(nr==1)
        if(any(strcmpi(style,{'ellip' 'cheby1'}))); pr=ripple;
        elseif(strcmpi(style,'cheby2')); sr=ripple;
        else sr=ripple; auto=1; 
        end
    elseif(nr==2)
        pr=ripple(1);
        sr=ripple(2);
        if(~strcmpi(style,'ellip')); auto=1; end
    else
        error('too many ripple parameters')
    end
end

% build and impliment filter for each group
for i=1:ng
    % sort and normalize corners
    pc=sort(corners(:)/nyq(i),1).';
    nc=length(pc);

    % check corners
    if(any(pc>=1)); error('corner(s) at/above nyquist');
    elseif(any(pc<=0)); error('corner(s) at/below 0'); end

    % dangerously close warnings
    if(any(pc>0.85))
        warning('SAClab:NyqClose','Corner close to Nyquist')
    end
    if(any(pc<0.001))
        warning('SAClab:ZeroClose','Extreme low frequency corner')
    end

    % work on corners by type
    if(any(strcmpi(type,{'low' 'high'})))
        % check number of corners
        if(~any(nc==[1 2]))
            error('incorrect number of corners defined')
        end
    
        % setup corners
        if(strcmpi(type,'low'))
            if(nc==1)
                if(~strcmpi(style,'cheby2'))
                    sc=pc*(1+tranbw);
                    if(sc>1); sc=0.9999; end
                else
                    sc=pc;
                    pc=pc*(1/(1+tranbw));
                end
            else
                sc=pc(2);
                pc=pc(1);
                auto=1;
            end
        else % high
            if(nc==1)
                if(~strcmpi(style,'cheby2'))
                    sc=pc-(1-pc)*tranbw;
                    if(sc<0); sc=0.0001; end
                else
                    sc=pc;
                    pc=sc+(1-sc)/(1/tranbw+1);
                end
            else
                sc=pc(1);
                pc=pc(2);
                auto=1;
            end
        end
    else % stop/notch or bandpass
        % check number of corners
        if(~any(nc==[2 4]))
            error('incorrect number of corners defined')
        end
    
        % setup corners
        if(nc==2)
            w=pc(2)-pc(1);
            if(~strcmpi(style,'cheby2'))
                sc(1)=pc(1)-tranbw*w;
                sc(2)=pc(2)+tranbw*w;
                if(sc(1)<0); sc(1)=0.0001; end
                if(sc(2)>1); sc(2)=0.9999; end
            else
                sc=pc;
                pc(1)=pc(1)+(1/(1/tranbw+2))*w;
                pc(2)=pc(2)-(1/(1/tranbw+2))*w;
            end
        else
            sc=[pc(1) pc(4)];
            pc=[pc(2) pc(3)];
            auto=1;
        end
    
        % swap if stop/notch
        if(strcmpi(type,{'stop'}))
            [pc,sc]=swap(sc,pc);
        end
    end
    
    % find filter order if needed
    if(auto==1)
        % determine by type
        if(strcmpi(style,'butter'))
            [order,pc]=buttord(pc,sc,pr,sr);
        elseif(strcmpi(style,'cheby1'))
            [order,pc]=cheb1ord(pc,sc,pr,sr);
        elseif(strcmpi(style,'cheby2'))
            [order,sc]=cheb2ord(pc,sc,pr,sr);
        else
            [order,pc]=ellipord(pc,sc,pr,sr);
        end
        disp(['NYQUIST: ' num2str(nyq(i)) ...
            'Hz  -->  AUTO FILTER ORDER: ' num2str(order)])
    end

    % make the filter
    if(strcmpi(style,'butter'))
        [z,p,k]=butter(order,pc);
    elseif(strcmpi(style,'cheby1'))
        [z,p,k]=cheby1(order,pr,pc);
    elseif(strcmpi(style,'cheby2'))
        [z,p,k]=cheby2(order,sr,sc);
    else
        [z,p,k]=ellip(order,pr,sr,pc);
    end

    % put into a filter structure
    [sos,g]=zp2sos(z,p,k);
    fs(i)=dfilt.df2tsos(sos,g);

    % implement
    if(passes==1)
        data(gi{i})=impfilt(data(gi{i}),fs(i),mirror);
    elseif(passes==2)
        data(gi{i})=impfilt(data(gi{i}),fs(i),mirror);
        data(gi{i})=imprevfilt(data(gi{i}),fs(i),mirror);
    elseif(passes==3)
        data(gi{i})=imprevfilt(data(gi{i}),fs(i),mirror);
    elseif(passes==4)
        data(gi{i})=imprevfilt(data(gi{i}),fs(i),mirror);
        data(gi{i})=impfilt(data(gi{i}),fs(i),mirror);
    end
end

end

function [data]=imprevfilt(data,fs,mirror)
%IMPREVFILT   Implements reverse pass filter
%   Implements a filter design on SAClab data records and makes appropriate
%   header updates.  Takes a mirror option which does pseudo-IC to limit 
%   edge effects.  Works with multiple records.

% combine records
[recs,ind,store,npts]=combo(data);
recs=recs(end:-1:1,:);

% mirror logical
if(mirror)
    % prepend a mirror-flip of the series to limit edge effects
    recs=filter(fs,[2*recs(ones(end-1,1),:)-recs(end:-1:2,:); recs]);
    recs=recs(fix(end/2)+1:end,:);
else
    % straight forward filter
    recs=filter(fs,recs);
end

% distribute records back
data=distro(data,recs(end:-1:1,:),ind,store,npts);

end

function [data]=impfilt(data,fs,mirror)
%IMPFILT   Implements filter
%   Implements a filter design on SAClab data records and makes appropriate
%   header updates.  Takes a mirror option which does pseudo-IC to limit 
%   edge effects.  Works with multiple records.

% combine records
[recs,ind,store,npts]=combo(data);

% mirror logical
if(mirror)
    % prepend a mirror-flip of the series to limit edge effects
    recs=filter(fs,[2*recs(ones(end-1,1),:)-recs(end:-1:2,:); recs]);
    recs=recs(fix(end/2)+1:end,:);
else
    % straight forward filter
    recs=filter(fs,recs);
end

% distribute records back
data=distro(data,recs,ind,store,npts);

end

