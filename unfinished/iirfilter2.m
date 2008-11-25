function [data,fs,nyq]=iirfilter2(data,type,style,varargin)
%IIRFILTER    Apply an IIR filter to SAClab data records
%
%    WARNING:  This code is in an unstable state.  Inputs, outputs, and
%     behavior will change with future versions.
%
%    Description: Utilizing the set of supplied parameters, IIR filters are
%     built to the desired specs and implemented on the data records.
%
%     Includes ability to automatically determine best fitting order
%     (number of poles) for the filter when the specs constain other
%     parameters like transition bandwidth, attenuation, etc.  Also by
%     default, filter implementation is done on a mirror-flipped expansion
%     of records to reduce edge-effects.
%
%
%    ORDER-IMPORTANT ARGUMENTS:
%
%     [DATA,FS,NYQ]=IIRFILTER(DATA,TYPE,STYLE,...)
%
%     DATA:   SAClab data structure
%
%     TYPE:   'high', 'low', 'notch'/'stop' or 'bandpass'
%
%     STYLE:  'butter', 'cheby1', 'cheby2' or 'ellip'
%
%             BUTTER:  Filter is monotonic in passband and stopband.
%             Any specs beyond passband corners and the filter order
%             require automatic filter order calculation to attempt to
%             implement them.  Order calculation will fit stopband corners
%             exactly and will always have a transition bandwidth less than
%             or equal to that designed.  This may lead to a wider passband
%             than desired.
%
%             CHEBY1: Filter is equiripple in the passband and monotonic in
%             the stopband.  Any specs beyond passband corners, maximum
%             passband attenuation, and the filter order requires automatic
%             filter order calculation.  Order calculation will fit 
%             passband corners exactly and will always have a transition 
%             bandwidth less than or equal to that designed.  This may lead 
%             to inaccurate stopband corner positioning.
%
%             CHEBY2: Filter is monotonic in the passband and equiripple in
%             the stopband.  Any specs beyond stopband corners, minimum
%             stopband attenuation, and the filter order requires automatic
%             filter order calculation.  Order calculation will fit 
%             stopband corners exactly and will always have a transition 
%             bandwidth less than or equal to that designed.  This may lead 
%             to a wider passband than desired.
%
%             ELLIP: Filter is equiripple in passband and stopband.
%             Automatic filter order determination is only required for
%             specs defining the stopband corner positioning.  Order 
%             calculation will fit passband corners exactly and will always
%             have a transition bandwidth less than or equal to that 
%             designed.  This may lead to inaccurate stopband corner 
%             positioning.
% 
%     FS:     Optional output. Filter structure.  Useful for visualizing 
%             and checking the filter with the fvtool.  Remember that the
%             filter(s) are tied to the sampling frequency of the records, 
%             so datasets with multiple samplerates will produce a 
%             structure with multiple filters!
%
%     NYQ:    Optional output. Nyquist frequencies associated with the 
%             filter(s) in the filter structure.
%
%
%    AT LEAST ONE OF THE FOLLOWING IS REQUIRED (SEE STYLE COMMENTS ABOVE 
%    FOR DETAILS - USING BOTH IS ALLOWED BUT SPECIFYING NPOLES IS PREFERRED
%    OVER THIS):
%
%     IIRFILTER(...,'PCORNER',CORNERS,...)
%     IIRFILTER(...,'SCORNER',CORNERS,...)
%
%     PCORNER: Array of pass corner(s).  Should be 1 element for low or
%              high pass filters and 2 elements for bandpass or notch
%              filters.  
%
%     PCORNERS: Array of pass corner(s) and optionally stop corner(s).  For
%              notch type filtering the pass corners are optional instead.  
%              If stop corner(s) are given, then the order parameter is
%              ignored and automatic filter order determination is done
%              utilizing those corners and the ripple parameters (if
%              given).  If no stop corner(s) are given and automatic order 
%              calculation is necessary (for instance if the ripple 
%              parameters are set), default stop corner(s) are created.  
%              These are set at 1/3 the passband width from the passband 
%              corner(s).  If this puts the corners outside the bandwidth 
%              of the signal (0Hz to nyquist) the corner is adjusted to be
%              near the nyquist/0Hz mark.
% 
%              Note that the automatic filter order determination will not
%              honor both the passband and stopband corners exactly as it
%              finds the best fitting order (remember that order is a
%              positive integer and so limits the design) that fits the 
%              corners given.  Usually the stopband corners are honored
%              over the passband corners, so if you want the pass corners 
%              honored it is best to just set the order and not give the 
%              stopband corners (although if you want to set the 
%              attenuation parameter to be non-default and don't find the
%              default corners satisfying you will have to set them).
%              There are checks in place that will warn when the designed
%              filter has corners that significantly differ from the ones
%              given.
%
%              Passband corners set the frequency at which the attenuation
%              hits 3dB, while stopband corners set the 30dB position.  The
%              dB levels that correspond to the corners can be adjusted
%              utilizing the attenuation parameter.
%
%              Lowpass or highpass filtering requires 1 or 2 inputs while
%              bandpass or notch filters require 2 or 4 inputs for the 
%              corners option.
%
%              ***  UNITS ARE IN Hz!!!  ***
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
%
%    OPTIONAL INPUTS:
%
%     IIRFILTER(...,'NPOLES',NUMBER_OF_POLES,...)
%     IIRFILTER(...,'PASSES',PASSES_OPTION,...)
%     IIRFILTER(...,'PATTEN',PASSBAND_ATTENUATION,...)
%     IIRFILTER(...,'SATTEN',STOPBAND_ATTENUATION,...)
%     IIRFILTER(...,'MIRROR',TRUE_OR_FALSE,...)
%     IIRFILTER(...,'TRANBW',FRACTION_OF_PASSBAND,...)
%
%     NPOLES: Optional but SUGGESTED.  If undefined or 0, automatic filter
%             order calculations are made (utilizes the supplied/default
%             corners and attenuation parameters to define a proper order).
%             Filter order is analogous to the sharpness/steepness/slope of
%             the transition from the passband to the stopband.  Higher 
%             orders provide better frequency resolution at the cost time 
%             resolution (aka ringing).  Note that automatic order
%             determination MAY ALTER CORNER POSITIONING SIGNIFICANTLY (a 
%             warning is issued if so).
%
%     PASSES: Optional.  Defaults to 1.  Accepts 1,2,3 or 4.  Option 1 is a
%             simple forward filter (single pass).  Option 2 filters 
%             forwards then backwards to perform zero-phase filtering.
%             Option 3 is for backwards filtering only and option 4 forward
%             filters after backwards filtering - aka a 2-pass in reverse 
%             order).  Two pass filtering is convenient for zero-phase 
%             filtering as single pass IIR filtering causes phase 
%             dispersion.  However, backwards filtering makes the signal 
%             acausal and should not be done in conjunction with onset 
%             picking.
%
%     Atten:  Optional.  1 or 2 element array of ripple/attenuation values
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
%             ***  UNITS IN dB!!!  ***
%
%    Usage: [data,fs,nyq]=iirfilter(data,type,style,corners,...
%                                        'order',order,...
%                                        'passes',passes,...
%                                        'atten',atten,...
%                                        'tranbw',tranbw,...
%                                        'mirror',mirror,...
%
%
%    Examples:
%      4th order lowpass butter filter with a passband corner of 10s
%     [data]=iirfilter(data,'low','butter',1/10,4)
%
%      bandpass elliptic filter with corners of 15s and 12s using
%      default stop corners/attenuation to determine filter order
%     [data]=iirfilter(data,'bandpass','ellip',[1/15 1/12])
%
%      2-pass 8th-order bandpass elliptic filter with 0.1dB passband ripple
%     [data]=iirfilter(data,'bandpass','ellip',[1/20 1/15],8,2,0.1)
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
%           unevenly spaced or xyz files will likely not be meaningufull
%           though!
%
%    See also: intrpol8, syncsr, deci, stretch, dft, idft

% check number of inputs
if(nargin<5); error('SAClab:iirfilter:badInput','not enough inputs'); end

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

% defaults
pc=[];      % pass corners (null)
sc=[];      % stop corners (null)
order=0;    % filter order (null)
auto=0;     % automatic order determination flag (off)
tranbw=1/3; % transition band / pass band width ratio
mirror=1;   % append mirror-flip flag (on)
pa=3;       % passband attenuation (in dB)
sa=30;      % stopband attenuation (in dB)
passes=1;   % default pass option (forward pass only)

% notch to stop
type=lower(type);
if(strcmp(type,'notch')); type='stop'; end

% check filter type and style
style=lower(style);
if(~any(strcmp(style,{'butter' 'cheby1' 'cheby2' 'ellip'})) || ...
   ~any(strcmp(type,{'low' 'high' 'stop' 'bandpass'})))
        error('SAClab:iirfilter:badInput',...
            'unknown filter: %s,%s',type,style)
end

% parse optional filter parameters
if(~mod(nargin,2))
    error('SAClab:iirfilter:unpairedInput',...
        'unpaired optional filter parameter'); 
end
for i=1:2:nargin-4
    switch lower(varargin{i})
        case 'pcorner'
            pc=sort(varargin{i+1}(:)).';
            if(any(pc<=0))
                error('SAClab:iirfilter:badInput',...
                    'CORNERS must be >0 (in Hz)')
            end
        case 'scorner'
            sc=sort(varargin{i+1}(:)).';
            if(any(sc<=0))
                error('SAClab:iirfilter:badInput',...
                    'CORNERS must be >0 (in Hz)')
            end
        case 'npoles'
            order=varargin{i+1};
            if(isempty(order))
                auto=1;
            elseif(~isscalar(order) || fix(order)~=order || order<0)
                error('SAClab:iirfilter:badInput',...
                    'NPOLES must be integer >=0')
            end
        case 'passes'
            passes=varargin{i+1};
            if(~isscalar(passes) || ~any(passes==[1 2 3 4]))
                error('SAClab:iirfilter:badInput',...
                    'unknown PASSES option')
            end
        case 'patten'
            pa=varargin{i+1};
            if(isempty(pa) || ~isscalar(pa) || pa<=0)
                error('SAClab:iirfilter:badInput',...
                    'ATTEN must be a scalar >0 (in dB)')
            end
        case 'satten'
            sa=varargin{i+1};
            if(isempty(sa) || ~isscalar(sa) || sa<=0)
                error('SAClab:iirfilter:badInput',...
                    'ATTEN must be a scalar >0 (in dB)')
            end
        case 'mirror'
            mirror=varargin{i+1};
        case 'tranbw'
            tranbw=varargin{i+1};
            if(~isscalar(tranbw) || tranbw<=0)
                error('SAClab:iirfilter:badInput',...
                    'TRANBW must be a scalar >0 (fraction of passband)')
            end
        otherwise
            error('SAClab:iirfilter:unknownOption',...
                'unknown option: %s',varargin{i})
    end
end

% check if order is 0 (auto order option)
if(order==0); auto=1; end

% check number of corners
npc=length(pc);
nsc=length(sc);
nc=npc+nsc;
if(any(strcmp(type,{'low' 'high'})))
    % check number of corners
    if(~any(npc==[0 1]) || ~any(nsc==[0 1]) || ~any(nc==[1 2]))
        error('SAClab:iirfilter:badInput',...
            'LOW or HIGH pass filter must have 1 or 2 corners defined')
    end
    
    % auto order if all corners given
    if(nc==2); auto=1; end
% stop/notch or bandpass
else
    % check number of corners
    if(~any(npc==[0 2]) || ~any(nsc==[0 2]) || ~any(nc==[2 4]))
        error('SAClab:iirfilter:badInput',...
            'BANDPASS or NOTCH filter must have 2 or 4 corners defined')
    end
    
    % auto order if all corners given
    if(nc==4); auto=1; end
end

% check attenuation setup (look if auto order option needed)
if(pa~=3)
    if(any(strcmp(style,{'butter' 'cheby2'}))); auto=1; end
end
if(sa~=30)
    if(any(strcmp(style,{'butter' 'cheby1'}))); auto=1; end
end

% make delta groups
delta=getheader(data,'delta');
gdt=unique(delta); ng=length(gdt); gi=cell(ng,1); nyq=1./(2*gdt);
for i=1:ng; gi{i}=find(delta==gdt(i)).'; end
if(ng>1)
    warning('SAClab:iirfilter:multipleSampleRates',...
        ['Dataset has multiple samplerates.  Filtering will proceed'...
        'by dividing the dataset into samplerate subsets and building'...
        'filters for each subset.  If samplerates differ significantly'...
        'the filters will likely differ significantly, which will lead'...
        'to unexpected results.'])
end

% build and impliment filter for each group
fs=cell(ng,1);
for i=1:ng
    % normalize corners
    normpc=pc/nyq(i);
    normsc=sc/nyq(i);

    % check corners
    if(any([normpc normsc]>=1))
        error('SAClab:iirfilter:badInput',...
            'corner(s) at/above nyquist')
    elseif(any([normpc normsc]<=0))
        error('SAClab:iirfilter:badInput',...
        'corner(s) at/below 0')
    end

    % frequency resolution warnings
    if(any([normpc normsc]>0.85))
        warning('SAClab:iirfilter:nyqClose',...
            'corner(s) close to Nyquist - consider upsampling!')
    end
    if(any([normpc normsc]<0.001))
        warning('SAClab:iirfilter:zeroClose',...
            'extreme low frequency corner(s) - consider decimating!')
    end

    % defaults corners
    if(strcmp(type,'low') && nc==1)
        % pass corner given, stop corner not
        if(npc==1)
            % find default stop corner
            normsc=normpc*(1+tranbw);
            
            % stop corner out of range
            % - reset corner to just inside 1
            % - sharper transition
            %if(normsc>=1); normsc=0.9999; end
        % stop corner given, pass corner not
        else
            % find default pass corner
            normpc=normsc*(1/(1+tranbw));
        end
    elseif(strcmp(type,'high') && nc==1)
        % pass corner given, stop corner not
        if(npc==1)
            % pass corner given, find default stop corner
            normsc=normpc-(1-normpc)*tranbw;
            
            % stop corner out of range
            % - reset corner to just inside 0
            % - sharper transition
            %if(normsc<=0); normsc=0.0001; end
        % stop corner given, pass corner not
        else
            % find default pass corner
            normpc=normsc+(1-normsc)/(1/tranbw+1);
        end
    elseif(strcmp(type,'bandpass') && nc==2)
        % pass corners given, stop corners not
        if(npc==2)
            % find default stop corners
            w=normpc(2)-normpc(1);
            normsc(1)=normpc(1)-tranbw*w;
            normsc(2)=normpc(2)+tranbw*w;
            
            % stop corners out of range
            % - reset corner to just inside 0 or 1
            % - sharper transition
            %if(normsc(1)<=0); normsc(1)=0.0001; end
            %if(normsc(2)>=1); normsc(2)=0.9999; end
        % stop corners given, pass corners not
        else
            % find default pass corners
            w=normsc(2)-normsc(1);
            normpc(1)=normsc(1)+(1/(1/tranbw+2))*w;
            normpc(2)=normsc(2)-(1/(1/tranbw+2))*w;
        end
    elseif(strcmp(type,'stop') && nc==2)
        % pass corners given, stop corners not
        if(npc==2)
            % find default stop corners
            w=normpc(2)-normpc(1);
            normsc(1)=normpc(1)+(1/(1/tranbw+2))*w;
            normsc(2)=normpc(2)-(1/(1/tranbw+2))*w;
        % stop corners given, pass corners not
        else
            % find default pass corners
            w=normsc(2)-normsc(1);
            normpc(1)=normsc(1)-tranbw*w;
            normpc(2)=normsc(2)+tranbw*w;
            
            % pass corners out of range
            % - reset corner to just inside 0 or 1
            % - sharper transition
            %if(normpc(1)<=0); normpc(1)=0.0001; end
            %if(normpc(2)>=1); normpc(2)=0.9999; end
        end
    end
    
    % find filter order if needed
    if(auto==1)
        % save corners (for check)
        normpc2=normpc;
        normsc2=normsc;
        
        % determine by type
        if(strcmp(style,'butter'))
            [order,normpc]=buttord(normpc,normsc,pa,sa);
        elseif(strcmp(style,'cheby1'))
            [order,normpc]=cheb1ord(normpc,normsc,pa,sa);
        elseif(strcmp(style,'cheby2'))
            [order,normsc]=cheb2ord(normpc,normsc,pa,sa);
        else
            [order,normpc]=ellipord(normpc,normsc,pa,sa);
        end
        
        % check that corners have not been significantly altered
        % - defining 'significant' as a disagreement of 0.01 or more
        if(any(abs([normpc2-normpc normsc2-normsc])>0.01))
            % you dun fcked up
            warning('SAClab:iirfilter:cornersNotPreserved',...
                ['\n\n!!!!! INACCURATE FILTER DESIGN !!!!!\n\n'...
                'Corners from automatic order determination\n'...
                'significantly differ from the desired corners.\n'...
                '\nPOSSIBLE REASONS:\n'...
                '1. You have designed a Butter filter with a\n'...
                '   non-default passband attenuation.\n'...
                '2. You have designed a weak filter (defined\n'...
                '   as a filter with a wide transition band)\n'...
                '   that is best fit by a order 1 or 2 filter.\n'...
                '\nThe first reason is a false alarm, but the problem\n'...
                'could actually be a combination of the two, so take\n'...
                'caution.  If you need some corners strictly enforced,\n'...
                'you will need to alter your supplied corners\n'...
                'or bypass automatic order calculation!\n\n']);
        end
        
        % display filter order
        disp(['NYQUIST: ' num2str(nyq(i)) ...
            'Hz  -->  AUTO FILTER ORDER: ' num2str(order)])
    end

    % make the filter
    if(strcmp(style,'butter'))
        [z,p,k]=butter(order,normpc,type);
    elseif(strcmp(style,'cheby1'))
        [z,p,k]=cheby1(order,pa,normpc,type);
    elseif(strcmp(style,'cheby2'))
        [z,p,k]=cheby2(order,sa,normsc,type);
    else
        [z,p,k]=ellip(order,pa,sa,normpc,type);
    end
    
    % put into a filter structure
    [sos,g]=zp2sos(z,p,k);
    fs(i)=dfilt.df2tsos(sos,g);
    
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
