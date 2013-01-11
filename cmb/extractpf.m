function [pf]=extractpf(pf,type,in1,in2)
%EXTRACTPF    Extracts profiles using the specified constraints
%
%    Usage:    pf=extractpf(pf,'freq',freq)
%              pf=extractpf(pf,'freq',freqrange)
%              pf=extractpf(pf,'ev',evdir)
%              pf=extractpf(pf,'st',stn)
%              pf=extractpf(pf,'st',stn1,stn2)
%
%    Description:
%     PF=EXTRACTPF(PF,'FREQ',FREQ) returns all profiles that have a
%     frequency range over the specified frequency FREQ.  FREQ is in Hz.
%
%     PF=EXTRACTPF(PF,'FREQ',FREQRANGE) returns all profiles that have a
%     frequency range within the specified frequency range FREQRANGE.
%     FREQRANGE is in Hz and must be [MINHZ MAXHZ].
%
%     PF=EXTRACTPF(PF,'EV',EVDIR) returns all profiles that originated from
%     waveforms in directory EVDIR.  EVDIR is just the waveform directory
%     name NOT the full path.
%
%     PF=EXTRACTPF(PF,'ST',STN) returns all profiles with the station
%     identified by STN.  This might be useful for inspecting the results
%     related to a particular station.  STN is either a string or a
%     cell-string identifying the station using one of the formats below:
%      (1) 'STATION'
%      (2) {'NETWORK' 'STATION'}
%      (3) {'NETWORK' 'STATION' 'STREAM'}
%      (4) {'NETWORK' 'STATION' 'STREAM' 'COMPONENT'}
%
%     PF=EXTRACTPF(PF,'ST',STN1,STN2) enables extracting profiles between a
%     specific pair of stations.
%
%    Notes:
%
%    Examples:
%     % Extract all profiles between CCM (Cathedral Cave, MO) & HRV
%     % (Harvard) -- please note that the input order is not important:
%     newpf=extractpf(pf,'st','CCM','HRV');
%
%     % Pull profiles for station CMB of the Berkeley network:
%     newpf=extractpf(pf,'st',{'BK' 'CMB'});
%
%    See also: SLOWDECAYPAIRS, PLOT_CMB_MEASUREMENTS, PLOT_CMB_PDF

%     Version History:
%        Feb.  2, 2011 - initial version
%        Mar. 30, 2011 - minor doc update
%        Mar.  2, 2012 - major rewrite, rename to extractpf
%        Oct. 11, 2012 - fix examples, drop corrections field requirement
%        Oct. 15, 2012 - st option works now, added more flexibility to
%                        option specification
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 15, 2012 at 13:35 GMT

% todo:

% check nargin
error(nargchk(3,4,nargin));

% check profile struct
reqfields={'gcdist','azwidth','slow','slowerr','decay','decayerr',...
    'cslow','cslowerr','cdecay','cdecayerr','cluster','kname','st','ev',...
    'delaz','synthetics','earthmodel','corrcoef','freq',...
    'phase','runname','dirname','time'};
if(~isstruct(pf) || any(~isfield(pf,reqfields)))
    error('seizmo:extractpf:badInput',...
        ['PF must be a struct with the fields:\n' ...
        sprintf('''%s'' ',reqfields{:}) '!']);
end

% check type
valid={'freq' 'fr' 'f' 'e' 'ev' 'eq' 's' 'st' 'stn' 'kname' 'kn' 'k' 'n'};
if(~ischar(type) || ~any(strcmpi(type,valid)))
    error('seizmo:extractpf:badInput',...
        'TYPE must be ''FREQ'' ''EV'' or ''ST''!');
end

% process by type
switch lower(type)
    case {'freq' 'fr' 'f'}
        switch numel(in1)
            case 1
                if(~isnumeric(in1) || ~isreal(in1) || in1<0)
                    error('seizmo:extractpf:badInput',...
                        'FREQ must be a positive real number!');
                end
                f=reshape([pf.freq],2,[]).';
                pf=pf(f(:,1)<=in1 & f(:,2)>=in1);
            case 2
                if(~isnumeric(in1) || ~isreal(in1) || any(in1<0) ...
                        || in1(1)>in1(2))
                    error('seizmo:extractpf:badInput',...
                        'FREQRANGE must be [MINHZ MAXHZ]!');
                end
                f=reshape([pf.freq],2,[]).';
                pf=pf(f(:,1)>=in1(1) & f(:,2)<=in1(2));
            otherwise
                error('seizmo:extractpf:badInput',...
                    'FREQ input invalid!');
        end
    case {'e' 'ev' 'eq'}
        if(~ischar(in1) || ndims(in1)>2 || size(in1,1)~=1)
            error('seizmo:extractpf:badInput',...
                'EVDIR must be a string!');
        end
        names={pf.dirname}';
        for i=1:numel(names)
            [a,b,c]=fileparts(names{i});
            names{i}=[b c];
        end
        pf=pf(ismember(names,in1));
    case {'st' 's' 'stn' 'kname' 'kn' 'k' 'n'}
        % check stn info
        if(isstring(in1)); in1=cellstr(in1); end
        if(~iscellstr(in1) || numel(in1)>4)
            error('seizmo:extractpf:badInput',...
                ['STN1 must be one of the following:\n' ...
                '''NAME'',{''NET'' ''NAME''},{''NET'' ' ...
                '''NAME'' ''STREAM''}']);
        end
        if(nargin==4)
            if(isstring(in2)); in2=cellstr(in2); end
            if(~iscellstr(in2) || numel(in2)>4)
                error('seizmo:extractpf:badInput',...
                    ['STN2 must be one of the following:\n' ...
                    '''NAME'',{''NET'' ''NAME''},{''NET'' ' ...
                    '''NAME'' ''STREAM''}']);
            end
        end
        
        % extract kname
        kname={pf.kname};
        nstns=cellfun('size',kname,1);
        
        % remove multistation profiles
        bad=nstns~=2;
        pf(bad)=[];
        kname(bad)=[];
        
        % quit if none
        if(isempty(pf))
            warning('seizmo:extractpf:no2StationFound',...
                'No 2-station profiles given!');
            return;
        end
        
        % combine kname info (each pair is a new page)
        kname=upper(reshape([kname{:}],2,4,[]));
        npf=size(kname,3);
        
        % find matches for first station
        match11=false(npf,1); match12=false(npf,1);
        in1=upper(in1); n1=numel(in1);
        switch n1
            case 1 % name
                match11(:)=strcmpi(in1,kname(1,2,:));
                match12(:)=strcmpi(in1,kname(2,2,:));
            case 2 % net name
                for i=1:npf
                    match11(i)=isequal(in1,kname(1,1:2,i));
                    match12(i)=isequal(in1,kname(2,1:2,i));
                end
            case 3 % net name hole
                for i=1:npf
                    match11(i)=isequal(in1,kname(1,1:3,i));
                    match12(i)=isequal(in1,kname(2,1:3,i));
                end
            case 4 % net name hole cmp
                for i=1:npf
                    match11(i)=isequal(in1,kname(1,1:4,i));
                    match12(i)=isequal(in1,kname(2,1:4,i));
                end
        end
        
        % find matches for station 2
        if(nargin==4)
            match21=false(npf,1); match22=false(npf,1);
            in2=upper(in2); n2=numel(in2);
            switch n2
                case 1 % name
                    match21(:)=strcmpi(in2,kname(1,2,:));
                    match22(:)=strcmpi(in2,kname(2,2,:));
                case 2 % net name
                    for i=1:npf
                        match21(i)=isequal(in2,kname(1,1:2,i));
                        match22(i)=isequal(in2,kname(2,1:2,i));
                    end
                case 3 % net name hole
                    for i=1:npf
                        match21(i)=isequal(in2,kname(1,1:3,i));
                        match22(i)=isequal(in2,kname(2,1:3,i));
                    end
                case 4 % net name hole cmp
                    for i=1:npf
                        match21(i)=isequal(in2,kname(1,1:4,i));
                        match22(i)=isequal(in2,kname(2,1:4,i));
                    end
            end
        end
        if(nargin==4)
            % find pair matches
            pf=pf((match11 & match22) | (match12 & match21));
        else
            % single station matches
            pf=pf(match11 | match12);
        end
end

end
