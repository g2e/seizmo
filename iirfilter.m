function [data,fs,nyq]=iirfilter(data,type,style,corners,order,passes,ripple)
%IIRFILTER    Apply an IIR filter to SEIZMO data records
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
%             and ripple/attenuation parameters to define a proper order).  
%             Filter order is analogous to the sharpness/steepness/slope of
%             the transition from the passband to the stopband.  Higher 
%             orders provide better frequency resolution at the cost time 
%             resolution (aka ringing).  Note that automatic order
%             determination HONORS THE STOPBAND CORNERS and may alter the
%             passband corners significantly to do this (a warning is
%             issued if so).
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
%     fs:     Optional output. Cell array of filter structure(s).  Useful
%             for visualizing the filter with the fvtool.  Remember that
%             the filter is tied to the sampling frequency of the records!
%
%     nyq:    Optional output. Nyquist frequencies for filter structure(s).
%
%
%    Usage: [data,fs]=iirfilter(data,type,style,corners,order,passes,ripple)
%
%
%    Examples:
%      4th order lowpass butter filter with a passband corner of 10s
%     [data]=iirfilter(data,'low','butter',1/10,4)
%
%      bandpass elliptic filter with corners of 15s and 12s using
%      default stop corners and ripples to determine filter order
%     [data]=iirfilter(data,'bandpass','ellip',[1/15 1/12])
%
%      2 pass bandpass elliptic filter with stop corners 
%      defined forcing order calc. 8th order is ignored but the
%      supplied 3dB passband ripple is included in the order calc.
%     [data]=iirfilter(data,'bandpass','ellip',[1/22 1/20 1/15 1/14],8,2,3)
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
%           Edge effect limiting in filtering is done by 'mirror-flipping' 
%           the data record before filtering as a pseudo-initial
%           conditioning measure.  The cost is that of filtering a record
%           twice as long.
%
%           This program will filter on any filetype.  The output for 
%           unevenly spaced or xyz files will not be accurate though!
%
%    See also: interpolate, syncrates, squish, stretch, dft, idft

% check number of inputs
error(nargchk(2,7,nargin))

% check data structure
error(seizmocheck(data,'dep'))

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% turn off header checking
oldcheckheaderstate=get_checkheader_state;
set_checkheader_state(false);

% check filter
if(~any(strcmpi(style,{'butter' 'cheby1' 'cheby2' 'ellip'})) || ...
   ~any(strcmpi(type,{'low' 'high' 'notch' 'stop' 'bandpass'})))
        error('seizmo:iirfilter:badInput','unknown filter')
end

% built-in defaults
auto=0;     % automatic order determination flag (off)
tranbw=1/3; % transition band / pass band width ratio
mirror=1;   % append mirror-flip flag (on)
pr=3;       % passband attenuation (in dB)
sr=30;      % stopband attenuation (in dB)
dp=1;       % default pass option (forward pass only)

% make delta groups
delta=getheader(data,'delta');
gdt=unique(delta); ng=length(gdt); gi=cell(ng,1); nyq=1./(2*gdt);
for i=1:ng; gi{i}=find(delta==gdt(i)).'; end

% notch to stop
if(strcmpi(type,'notch')); type='stop'; end

% check order
if(nargin==4 || isempty(order) || order==0); auto=1; end

% check passes
if(nargin<6 || isempty(passes)); passes=dp;
elseif (~any(passes==[1 2 3 4]))
    error('seizmo:iirfilter:badInput','1 or 2 passes only')
end

% check ripple
if(nargin==7 && ~isempty(ripple))
    ripple=sort(ripple(:).');
    nr=length(ripple);

    % check for negative
    if(any(ripple<=0))
        error('seizmo:iirfilter:badInput',...
            'negative or zero ripple not allowed')
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
        error('seizmo:iirfilter:badInput',...
            'too many ripple parameters')
    end
end

% build and impliment filter for each group
fs=cell(ng,1);
for i=1:ng
    % sort and normalize corners
    pc=sort(corners(:)/nyq(i),1).';
    nc=length(pc);

    % check corners
    if(any(pc>=1))
        error('seizmo:iirfilter:badInput',...
            'corner(s) at/above nyquist')
    elseif(any(pc<=0))
        error('seizmo:iirfilter:badInput',...
        'corner(s) at/below 0')
    end

    % dangerously close warnings
    if(any(pc>0.85))
        warning('seizmo:iirfilter:nyqClose',...
            'Corner close to Nyquist')
    end
    if(any(pc<0.001))
        warning('seizmo:iirfilter:zeroClose',...
            'Extreme low frequency corner')
    end

    % work on corners by type
    if(any(strcmpi(type,{'low' 'high'})))
        % check number of corners
        if(~any(nc==[1 2]))
            error('seizmo:iirfilter:badInput',...
                'incorrect number of corners defined')
        end
    
        % setup corners
        if(strcmpi(type,'low'))
            if(nc==1)
                if(~strcmpi(style,'cheby2'))
                    sc=pc*(1+tranbw);
                    if(sc>=1); sc=0.9999; end
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
                    if(sc<=0); sc=0.0001; end
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
            error('seizmo:iirfilter:badInput',...
                'incorrect number of corners defined')
        end
    
        % setup corners
        if(nc==2)
            w=pc(2)-pc(1);
            if(~strcmpi(style,'cheby2'))
                sc(1)=pc(1)-tranbw*w;
                sc(2)=pc(2)+tranbw*w;
                if(sc(1)<=0); sc(1)=0.0001; end
                if(sc(2)>=1); sc(2)=0.9999; end
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
        % save pass corners
        pc2=pc;
        
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
        
        % check that passband range has not significantly altered
        % - defining 'significant' as a disagreement of 0.01 or more
        if(~strcmpi(style,'cheby2'))
            if(any(abs(pc2-pc)>0.01))
                % you dun fcked up
                warning('seizmo:iirfilter:cornersNotPreserved',...
                    ['\n\n!!!!! BAD FILTER DESIGN !!!!!\n\n'...
                    'Passband corners from automatic order determination\n'...
                    'significantly differ from the desired passband\n'...
                    'corners in order to satisfy the stopband corners.\n'...
                    'If you need the passband corners strictly enforced,\n'...
                    'you will need to alter the stopband corners, filter\n'...
                    'order and/or ripple/attenuation parameters!\n\n']);
            end
        end
        
        disp(['NYQUIST: ' num2str(nyq(i)) ...
            'Hz  -->  AUTO FILTER ORDER: ' num2str(order)])
    end

    % make the filter
    if(strcmpi(style,'butter'))
        [z,p,k]=butter(order,pc,type);
    elseif(strcmpi(style,'cheby1'))
        [z,p,k]=cheby1(order,pr,pc,type);
    elseif(strcmpi(style,'cheby2'))
        [z,p,k]=cheby2(order,sr,sc,type);
    else
        [z,p,k]=ellip(order,pr,sr,pc,type);
    end

    % put into a filter structure
    [sos,g]=zp2sos(z,p,k);
    fs{i}=dfilt.df2tsos(sos,g);
    
    % combine records
    [recs,idx1,ind,idx2,store,npts]=combinerecords(data(gi{i}));

    % implement
    if(passes==1)
        recs=impfilt(recs,fs{i},mirror,true);
    elseif(passes==2)
        recs=impfilt(impfilt(recs,fs{i},mirror,true),fs{i},mirror,false);
    elseif(passes==3)
        recs=impfilt(recs,fs{i},mirror,false);
    elseif(passes==4)
        recs=impfilt(impfilt(recs,fs{i},mirror,false),fs{i},mirror,true);
    end
    
    % distribute records back
    data(gi{i})=distributerecords(data(gi{i}),recs,idx1,[],[],store,npts);
end

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end

function [recs]=impfilt(recs,fs,mirror,forward)
%IMPREVFILT   Implements filter
%   Implements a filter design on SEIZMO data records and makes appropriate
%   header updates.  Takes a mirror option which does pseudo-IC to limit 
%   edge effects.  Works with multiple records.

% reverse
if(~forward); recs=recs(end:-1:1,:); end

% mirror logical
if(mirror)
    % prepend a mirror-flip of the series to limit edge effects
    recs=filter(fs,[2*recs(ones(end-1,1),:)-recs(end:-1:2,:); recs]);
    recs=recs(ceil(end/2):end,:);
else
    % straight forward filter
    recs=filter(fs,recs);
end

% reverse
if(~forward); recs=recs(end:-1:1,:); end

end
