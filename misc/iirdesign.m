function [fo,passes,mirror]=iirdesign(nyq,type,style,varargin)
%IIRDESIGN    Designs an iir filter with the given constraints
%
%    Usage:    [fo,passes,mirror]=...
%                    iirdesign(nyq,type,style,'option1',value1,...)
%              [fo,passes,mirror]=...
%                    iirdesign(nyq,fo,nyq2,'option1',value1,...)
%
%    Description:
%     [FO,PASSES,MIRROR]=IIRDESIGN(NYQ,TYPE,STYLE,'OPTION1',VALUE1,...)
%     designs an iir (infinite impulse response) filter using the design
%     specs specified.  The filter is then normalized for each Nyquist
%     frequency given in NYQ.  The filter(s) are output in cell array FO
%     which contains the corresponding dfilt.df2tsos filter object for each
%     Nyquist frequency.  The design process includes the ability to auto-
%     determine the best fitting order (aka number of poles) when the
%     filter is over designed.  Usually this is when the specs constain
%     parameters like transition bandwidth, stopband corners, attenuation,
%     etc.  See the options below for details/defaults.  Use the function
%     FVTOOL to verify/visualize the filter(s).
%
%      Valid TYPE strings (NO DEFAULT):
%       'LO' - Low-pass filter
%       'HI' - High-pass filter
%       'BP' - Band-pass filter
%       'BS' - Band-stop filter
%
%      Valid STYLE strings (NO DEFAULT):
%       'BUTTER':  Filter is monotonic in passband and stopband.
%             Any specs beyond passband corner(s) and the filter order
%             triggers automatic filter order calculation to attempt to
%             implement them.  Order calculation will fit passband corners
%             exactly and will always have a transition bandwidth less than
%             or equal to that designed.  This may lead to inaccurate
%             stopband corner positioning.
%
%       'CHEBY1':  Filter is equiripple in the passband and monotonic in
%             the stopband.  Any specs beyond passband corners, maximum
%             passband attenuation, and the filter order triggers automatic
%             filter order calculation.  Order calculation will fit
%             passband corners exactly and will always have a transition
%             bandwidth less than or equal to that designed.  This may lead
%             to inaccurate stopband corner positioning.
%
%       'CHEBY2':  Filter is monotonic in the passband and equiripple in
%             the stopband.  Any specs beyond stopband corners, minimum
%             stopband attenuation, and the filter order triggers automatic
%             filter order calculation.  Order calculation will fit 
%             stopband corners exactly and will always have a transition 
%             bandwidth less than or equal to that designed.  This will 
%             give a wider passband than desired.
%
%       'ELLIP':  Filter is equiripple in passband and stopband.
%             Automatic filter order determination is only required for
%             specs defining the stopband corner positioning.  Order 
%             calculation will fit passband corners exactly and will always
%             have a transition bandwidth less than or equal to that 
%             designed.  This may lead to inaccurate stopband corner 
%             positioning.
%
%      AT LEAST ONE OF THE FOLLOWING OPTIONS IS REQUIRED:
%       IIRDESIGN(...,'PCORNER',CORNERS,...)
%       IIRDESIGN(...,'SCORNER',CORNERS,...)
%
%      THE REST ARE OPTIONAL:
%       IIRDESIGN(...,'NPOLES',NUMBER_OF_POLES,...)
%       IIRDESIGN(...,'PASSES',PASSES_OPTION,...)
%       IIRDESIGN(...,'PATTEN',PASSBAND_ATTENUATION,...)
%       IIRDESIGN(...,'SATTEN',STOPBAND_ATTENUATION,...)
%       IIRDESIGN(...,'MIRROR',TRUE_OR_FALSE,...)
%       IIRDESIGN(...,'TRANBW',FRACTION_OF_PASSBAND,...)
%
%      PCORNER: Array of passband corner(s). Must be 1 element for low/high
%       pass filters and 2 elements for bandpass/stop filters.  There is no
%       default value for PCORNER, but if SCORNER is defined then PCORNER
%       may be automatically defined using SCORNER and TRANBW.  PCORNER may
%       be set to an empty array if such behavior is desired.  Note that
%       PCORNER will be honored in auto-order determination for 'BUTTER',
%       'CHEBY1' & 'ELLIP' style filters.  'CHEBY2'-style filters only
%       honors the stopband corners (the passband will thus be wider than
%       requested).
%
%      SCORNER: Array of stopband corner(s). Must be 1 element for low/high
%       pass filters and 2 elements for bandpass/stop filters.  There is no
%       default value for SCORNER, but if PCORNER is defined then SCORNER
%       may be automatically defined using PCORNER and TRANBW.  SCORNER may
%       be set to an empty array if such behavior is desired.  Note that
%       SCORNER will be honored only for 'CHEBY2'-style filters (the
%       attenuation will be stronger than requested for the other styles).
%
%      NPOLES: Optional but SUGGESTED!  If undefined or 0, automatic filter
%       order calculations are made (utilizes the supplied/default corners
%       and attenuation parameters to define a suitable order).  Filter
%       order is analogous to the sharpness/steepness/slope of the
%       transition from the passband to the stopband.  Higher orders
%       provide better frequency resolution at the cost time resolution
%       (aka ringing).  Note that automatic order determination MAY ALTER
%       CORNER POSITIONING SIGNIFICANTLY (a warning is issued).
%
%      PASSES: Default is 1.  Accepts 1, 2, -1 or -2.  Option 1 is a simple
%       forward filter (single pass).  Option 2 filters forwards then
%       backwards to perform zero-phase filtering.  Option -1 is for
%       backwards filtering only and option -2 forward filters after
%       backwards filtering - aka a 2-pass in reverse order).  Two-pass
%       filtering is convenient for zero-phase filtering as a single-pass
%       IIR filter causes phase/group dispersion.  However, backwards
%       filtering makes the resulting signal acausal and should not be done
%       in studies where the onset of a phase needs to be preserved.  Does
%       not affect the design, just the implementation (hence it is an
%       output argument).
%
%      PATTEN: Default is 3dB.  Indicates the maximum permissible passband
%       amplitude loss/ripple.
%
%      SATTEN: Default is 30dB.  Indicates the minimum permissible stopband
%       loss or the maximum ripple (ie. nothing in the stopband is allowed
%       above this value).
%
%      MIRROR: Default is TRUE.  Indicates whether the filter is
%       implemented on record(s) with a mirror-flip leading the record
%       (TRUE) or if filtering is done in the normal fashion (FALSE).
%       Attaching a mirror-flip to the beginning of a record prior to
%       filtering is a method of reducing edge-effects by essentially
%       removing the jump from zero at the start of the record.  The cost
%       of this operation is that of filtering a record approximately twice
%       as long.  Does not affect the design, just the implementation
%       (hence it is an output argument).
%
%      TRANBW: Default is 1/3.  Defines the width of the transition from
%       the stopband to the passband in terms of fraction of the passband
%       width.  Defining TRANBW will indicated the location of the
%       complimentary corner(s).  If both the passband and stopband
%       corner(s) are defined, then TRANBW is ignored.  Setting TRANBW
%       (even to the default value) will trigger automatic filter order
%       calculation.  To avoid this, set TRANBW to an empty array or do not
%       call the option.
%
%     [FO,PASSES,MIRROR]=...
%                    IIRDESIGN(NYQ,FO,NYQ2,'OPTION1',VALUE1,...)
%      allows for bypassing the design process (the filters are given in FO
%      and their corresponding Nyquist frequencies are in NYQ2).  This call
%      style will check that all the required Nyquist frequencies in NYQ
%      are found in NYQ2, rearranging FO to match NYQ.  Optional arguments
%      are parsed but no filter alteration/design is performed.  This is
%      mainly useful for SEIZMO function IIRFILTER.
%
%    Notes: For narrow band filters, PARTICULARLY AT LOW Hz, consider 
%           decimating the signal (using SQUISH or SYNCRATES) prior to
%           filtering to aid in the accuracy of the filter design.  This is
%           because the analog to digital transformation of the filter
%           design will lose accuracy if the features of the filter are
%           too narrow.  Basically try to keep the filter width above a few
%           percent of the total bandwidth of the signal.  
%
%           Remember, length of the trace does not affect filter design,
%           but the sampling/Nyquist frequency does.
%
%    Examples:
%     % 4th order lowpass butter filter with a passband corner of 10s:
%     fo=iirdesign(1,'lo','butter','c',1/10,'np',4)
%
%     % Bandpass elliptic filter with passband corners of 15s and 12s using
%     % default transition bandwidth & attenuation to get the filter order:
%     fo=iirdesign(1,'bp','ellip','c',[1/15 1/12])
%
%     % 2-pass 4th-order bandpass elliptic filter w/ 0.1dB passband ripple:
%     [fo,npass]=iirdesign(1,'bp','e','c',[1/20 1/15],'o',4,'p',2,'pr',0.1)
%
%     % Verify filters with the FVTOOL (Signal Processing Toolbox):
%     fo=iirdesign([10 20],'bp','e','c',[1/20 1/15],'o',4);
%     fvtool(fo{:}) % note that this will plot all the filters together
%
%    See also: IIRFILTER, FVTOOL, BUTTER, CHEBY1, CHEBY2, ELLIP,
%              BUTTORD2, CHEB1ORD, CHEB2ORD, ELLIPORD

%     Version History:
%        Feb.  2, 2010 - initial version
%        Mar. 25, 2010 - minor doc update
%        Apr. 25, 2010 - allow 'lp'/'hp' as filter types
%        May  20, 2010 - allow 'n' to specify poles, 'bu' for butter
%        Sep. 20, 2010 - alter 3,4 pass to -1,-2 for simplicity
%        Feb. 11, 2011 - mass nargchk fix
%        Feb.  4, 2012 - doc update, fixed formatting of a couple warnings
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  4, 2012 at 15:05 GMT

% todo:

% check number of inputs
error(nargchk(2,inf,nargin));

% defaults
default.ORDER=0;    % filter order (null)
default.TRANBW=1/3; % transition band / passband width ratio
default.PA=3;       % passband attenuation (in dB)
default.SA=30;      % stopband attenuation (in dB)
default.AUTO=0;     % automatic order determination flag (off)
default.MIRROR=1;   % append mirror-flip flag (on)
default.PASSES=1;   % default pass option (forward pass only)

% preallocate as empty
pc=[]; sc=[]; order=[]; tranbw=[]; auto=[];
mirror=[]; pa=[]; sa=[]; passes=[];

% check nyquist frequencies
if(isempty(nyq))
    error('seizmo:iirdesign:noNyquist',...
        'No Nyquist frequency given!');
elseif(~isreal(nyq) || any(nyq(:)<=0))
    error('seizmo:iirdesign:badNyquist',...
        'Nyquist frequencies are not positive real!');
end

% look out for dfilt.df2tsos object
dfo=false;
if(isa(type,'dfilt.df2tsos'))
    fo={type};
    dfo=true;
elseif(iscell(type))
    dfo=true;
    for i=1:numel(type)
        if(~isa(type{i},'dfilt.df2tsos'))
            dfo=false;
        end
    end
    if(dfo); fo=type; end
end

% dfilt object vs design
if(dfo)
    nyq2=style;
    % need to check/reorder fs to match nyq
    [ok,idx]=ismember(nyq,nyq2);
    if(~all(ok))
        error('seizmo:iirdesign:missingFilters',...
            ['Missing filter object(s) for some sample rates!' ...
            '\nAvailable Sample Rates: ' sprintf('%g ',2*nyq2) ...
            '\nMissing Sample Rates: ' sprintf('%g ',2*nyq(~ok))]);
    end
    fo=fo(idx);
else
    % check type
    if(~ischar(type))
        error('seizmo:iirdesign:badType',...
            'TYPE must be a string!');
    end
    switch lower(type)
        case {'low' 'lo' 'l' 'lp'}
            type='low';
        case {'high' 'hi' 'h' 'hp'}
            type='high';
        case {'bandpass' 'pass' 'bp'}
            type='bandpass';
        case {'bandstop' 'stop' 'bs' 'notch' 'n' ...
                'bandreject' 'reject' 'br'}
            type='stop';
        otherwise
            error('seizmo:iirdesign:badType',...
                'TYPE is unknown: %s',type);
    end

    % check style
    if(~ischar(style))
        error('seizmo:iirdesign:badStyle',...
            'STYLE must be a string!');
    end
    switch lower(style)
        case {'butter' 'butt' 'butterworth' 'bu' 'b'}
            style='butter';
        case {'cheby1' 'cheb1' 'chebychev1' 'c1'}
            style='cheby1';
        case {'cheby2' 'cheb2' 'chebychev2' 'c2'}
            style='cheby2';
        case {'ellip' 'el' 'elliptical' 'e'}
            style='ellip';
        otherwise
            error('seizmo:iirdesign:badStyle',...
                'STYLE is unknown: %s',style);
    end
end

% parse optional filter parameters
if(mod(nargin-3,2))
    error('seizmo:iirdesign:unpairedInput',...
        'Unpaired OPTION/VALUE!');
end
for i=1:2:nargin-3
    if(~ischar(varargin{i}))
        error('seizmo:iirdesign:badOption',...
            'Argument %d: OPTION must be a string!',i+3);
    end
    if(isempty(varargin{i+1})); continue; end
    switch lower(varargin{i})
        case {'auto' 'au'}
            if(~isscalar(varargin{i+1}) ...
                    || (~islogical(varargin{i+1}) ...
                    && ~any(varargin{i+1}==[0 1])))
                error('seizmo:iirdesign:badValue',...
                    'AUTO must be a scalar 0 (false) or 1 (true)!');
            end
            auto=varargin{i+1};
        case {'pcorners' 'pcorner' 'pc' 'corners' 'corner' 'c'}
            if(~any(numel(varargin{i+1})==[1 2]) ...
                    || ~isreal(varargin{i+1}) ...
                    || any(varargin{i+1}<=0))
                error('seizmo:iirdesign:badInput',...
                    'Passband CORNERS must be 1 to 2 reals >0 (in Hz)!');
            end
            pc=sort(varargin{i+1}(:)).';
        case {'scorners' 'scorner' 'sc'}
            if(~any(numel(varargin{i+1})==[1 2]) ...
                    || ~isreal(varargin{i+1}) ...
                    || any(varargin{i+1}<=0))
                error('seizmo:iirdesign:badInput',...
                    'Stopband CORNERS must be 1 to 2 reals >0 (in Hz)!');
            end
            sc=sort(varargin{i+1}(:)).';
        case {'npoles' 'order' 'np' 'n' 'o'}
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) ...
                    || fix(varargin{i+1})~=varargin{i+1} ...
                    || varargin{i+1}<0)
                error('seizmo:iirdesign:badInput',...
                    'ORDER must be an integer >=0 !');
            end
            order=varargin{i+1};
        case {'passes' 'pass' 'p'}
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) ...
                    || ~any(abs(varargin{i+1})==[1 2]))
                error('seizmo:iirdesign:badInput',...
                    'PASSES must be 1, 2, -1 or -2 !');
            end
            passes=varargin{i+1};
        case {'patten' 'pa' 'pripple' 'pr' 'ripple' 'r'}
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) ...
                    || varargin{i+1}<=0)
                error('seizmo:iirdesign:badInput',...
                    'Passband RIPPLE/ATTEN must be a scalar >0 (in dB)!');
            end
            pa=varargin{i+1};
        case {'satten' 'sa' 'sripple' 'sr' 'atten' 'at'}
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) ...
                    || varargin{i+1}<=0)
                error('seizmo:iirdesign:badInput',...
                    'Stopband RIPPLE/ATTEN must be a scalar >0 (in dB)!');
            end
            sa=varargin{i+1};
        case {'mirror' 'm'}
            if(~isscalar(varargin{i+1}) ...
                    || (~islogical(varargin{i+1}) ...
                    && ~any(varargin{i+1}==[0 1])))
                error('seizmo:iirdesign:badValue',...
                    'MIRROR must be a scalar 0 (false) or 1 (true)!');
            end
            mirror=varargin{i+1};
        case {'tranbw' 'tr' 't'}
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) ...
                    || varargin{i+1}<=0)
                error('seizmo:iirdesign:badInput',...
                    'TRANBW must be a scalar >0 (fraction of passband)!')
            end
            tranbw=varargin{i+1};
        otherwise
            error('seizmo:iirdesign:unknownOption',...
                'Unknown Option: %s',varargin{i})
    end
end

% check auto, mirror, passes, order
if(isempty(auto)); auto=default.AUTO; end
if(isempty(mirror)); mirror=default.MIRROR; end
if(isempty(passes)); passes=default.PASSES; end
if(isempty(order)); order=default.ORDER; end

% quick exit if dfilt object given
if(dfo); return; end

% check if order is 0 (auto order option)
if(~order); auto=1; end

% type specific checks (corner counting)
npc=numel(pc);
nsc=numel(sc);
nc=npc+nsc;
switch type
    case {'low' 'high'}
        % check number of corners
        if(~any(npc==[0 1]) || ~any(nsc==[0 1]) || ~any(nc==[1 2]))
            error('seizmo:iirdesign:notEnoughCorners',...
                ['Low/High pass filters may only have\n' ...
                '1 passband and/or 1 stopband corner defined!']);
        end
        
        % auto order if all corners given
        if(nc==2); auto=1; end
    case {'bandpass' 'stop'}
        % check number of corners
        if(~any(npc==[0 2]) || ~any(nsc==[0 2]) || ~any(nc==[2 4]))
            error('seizmo:iirdesign:notEnoughCorners',...
                ['Band pass/stop filter may only have\n' ...
                '2 passband and/or 2 stopband corners defined!']);
        end
        
        % auto order if all corners given
        if(nc==4); auto=1; end
end

% style specific checks (over designed filter)
switch style
    case 'butter'
        % 2 pc & order
        if(~isempty([pa sa sc tranbw])); auto=1; end
    case 'cheby1'
        % 2 pc & order & pa
        if(~isempty([sa sc tranbw])); auto=1; end
    case 'cheby2'
        % 2 sc & order & sa
        if(~isempty([pa pc tranbw])); auto=1; end
    case 'ellip'
        % 2 pc & order & pa & sa
        if(~isempty([sc tranbw])); auto=1; end
end

% fill in parameter gaps
[pc,sc,pa,sa]=getdesignargs(default,pc,sc,pa,sa,tranbw,nc,npc,type);

% build filters
nnyq=numel(nyq);
fo=cell(nnyq,1);
for i=1:nnyq
    % normalize corners
    normpc=pc/nyq(i);
    normsc=sc/nyq(i);
    
    % find filter order if needed
    if(auto)
        % check within (0,1)
        if(any([normpc normsc]>=1))
            error('seizmo:iirdesign:badCorner',...
                ['Corner(s) at/above Nyquist!\n' ...
                'Nyquist: %gHz\nPassband Corner(s): ' ...
                sprintf('%gHz ',normpc*nyq(i)) ...
                '\nStopband Corner(s): ' ...
                sprintf('%gHz ',normsc*nyq(i))],nyq(i));
        elseif(any([normpc normsc]<=0))
            error('seizmo:iirdesign:badCorner',...
                ['Corner(s) at/below 0Hz!' ...
                '\nPassband Corner(s): ' ...
                sprintf('%gHz ',normpc*nyq(i)) ...
                '\nStopband Corner(s): ' ...
                sprintf('%gHz ',normsc*nyq(i))]);
        end
        
        % save corners (for check)
        normpc2=normpc;
        normsc2=normsc;
        
        % determine by type
        switch style
            case 'butter'
                % buttord - transition band will partly be at/above pa
                % buttord2 - transition band will partly be at/below sa
                [order,normpc]=buttord2(normpc,normsc,pa,sa);
            case 'cheby1'
                % transition band will partly be at/below sa
                [order,normpc]=cheb1ord(normpc,normsc,pa,sa);
            case 'cheby2'
                % transition band will partly be at/above pa
                [order,normsc]=cheb2ord(normpc,normsc,pa,sa);
            case 'ellip'
                % transition band will partly be at/below sa
                [order,normpc]=ellipord(normpc,normsc,pa,sa);
        end
        
        % check that corners have not been altered
        if(~isequal(normpc,normpc2))
            % let the user know
            disp('PASSBAND CORNER(S) ADJUSTED!');
            for j=1:numel(normpc)
                fprintf('%gHz  -->  %gHz\n',...
                    normpc2(j)*nyq(i),normpc(j)*nyq(i));
            end
            disp('Verify the filter is acceptable with FVTOOL!');
        end
        if(~isequal(normsc,normsc2))
            % let the user know
            disp('STOPBAND CORNER(S) ADJUSTED!');
            for j=1:numel(normsc)
                fprintf('%gHz  -->  %gHz\n',...
                    normsc2(j)*nyq(i),normsc(j)*nyq(i));
            end
            disp('Verify the filter is acceptable with FVTOOL!');
        end
        
        % display filter order
        fprintf('NYQUIST: %gHz --> AUTO FILTER ORDER: %d\n',nyq(i),order);
    end

    % which corners to check
    corners=normpc; if(strcmp(type,'cheby2')); corners=normsc; end
    
    % check within (0,1)
    if(any(corners>=1))
        error('seizmo:iirdesign:badCorner',...
            ['Corner(s) at/above Nyquist!\n' ...
            'Nyquist: %gHz' ...
            '\nCorner(s): ' sprintf('%gHz ',corners*nyq(i))],nyq(i));
    elseif(any(corners<=0))
        error('seizmo:iirdesign:badCorner',...
            ['Corner(s) at/below 0Hz!' ...
            '\nCorner(s): ' sprintf('%gHz ',corners*nyq(i))]);
    end

    % frequency resolution warnings
    if(any(corners>0.85))
        warning('seizmo:iirdesign:nyqClose',...
            ['Corner(s) close to Nyquist - consider upsampling!' ...
            '\nNyquist: %gHz\nCorner(s): ' ...
            sprintf('%gHz ',corners*nyq(i))],nyq(i));
    end
    if(any(corners<0.001))
        warning('seizmo:iirdesign:zeroClose',...
            ['Extreme low frequency corner(s) - consider decimating!' ...
            '\nNyquist: %gHz\nCorner(s): ' ...
            sprintf('%gHz ',corners*nyq(i))],nyq(i));
    end
    
    % make the filter
    switch style
        case 'butter'
            [z,p,k]=butter(order,corners,type);
        case 'cheby1'
            [z,p,k]=cheby1(order,pa,corners,type);
        case 'cheby2'
            [z,p,k]=cheby2(order,sa,corners,type);
        case 'ellip'
            [z,p,k]=ellip(order,pa,sa,corners,type);
    end
    
    % put into a filter structure
    [sos,g]=zp2sos(z,p,k);
    fo{i}=dfilt.df2tsos(sos,g);
end

end

function [pc,sc,pa,sa]=getdesignargs(def,pc,sc,pa,sa,tranbw,nc,npc,type)
%GETDESIGNARGS    Returns parameters for filter design
%
% Basically:
% - grabs defaults for pa,sa,tranbw if needed
% - finds pc or sc given that the other is defined

% make sure we have attenuation/ripple parameters
if(isempty(pa)); pa=def.PA; end
if(isempty(sa)); sa=def.SA; end

% make sure transition bandwidth is set
if(isempty(tranbw)); tranbw=def.TRANBW; end

% defaults corners
if(strcmp(type,'low') && nc==1)
    % pass corner given, stop corner not
    if(npc)
        % find default stop corner
        sc=pc*(1+tranbw);
    % stop corner given, pass corner not
    else
        % find default pass corner
        pc=sc*(1/(1+tranbw));
    end
elseif(strcmp(type,'high') && nc==1)
    % pass corner given, stop corner not
    if(npc)
        % pass corner given, find default stop corner
        sc=pc-(1-pc)*tranbw;
    % stop corner given, pass corner not
    else
        % find default pass corner
        pc=sc+(1-sc)/(1/tranbw+1);
    end
elseif(strcmp(type,'bandpass') && nc==2)
    % pass corners given, stop corners not
    if(npc)
        % find default stop corners
        w=pc(2)-pc(1);
        sc(1)=pc(1)-tranbw*w;
        sc(2)=pc(2)+tranbw*w;
    % stop corners given, pass corners not
    else
        % find default pass corners
        w=sc(2)-sc(1);
        pc(1)=sc(1)+(1/(1/tranbw+2))*w;
        pc(2)=sc(2)-(1/(1/tranbw+2))*w;
    end
elseif(strcmp(type,'stop') && nc==2)
    % pass corners given, stop corners not
    if(npc)
        % find default stop corners
        w=pc(2)-pc(1);
        sc(1)=pc(1)+(1/(1/tranbw+2))*w;
        sc(2)=pc(2)-(1/(1/tranbw+2))*w;
    % stop corners given, pass corners not
    else
        % find default pass corners
        w=sc(2)-sc(1);
        pc(1)=sc(1)-tranbw*w;
        pc(2)=sc(2)+tranbw*w;
    end
end

end
