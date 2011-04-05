function [pf]=extractpair(pf,stn1,stn2)
%EXTRACTPAIR    Extracts profiles between the specified station pair
%
%    Usage:    pf=extractpair(pf,stn1,stn2)
%              pf=extractpair(pf,stn)
%
%    Description:
%     PF=EXTRACTPAIR(PF,STN1,STN2) enables extracting profiles between a
%     specific pair of stations from the core-diffracted profile struct PF.
%     STN1 & STN2 are either strings or cell-strings identifying the
%     stations using one of the formats below:
%      (1) 'STATION'
%      (2) {'NETWORK' 'STATION'}
%      (3) {'NETWORK' 'STATION' 'STREAM'}
%
%     PF=EXTRACTPAIR(PF,STN) returns all 2-station profiles with the
%     station identified by STN.  This might be useful for inspecting the
%     results related to a station.
%
%    Notes:
%
%    Examples:
%     % Extract all profiles between CCM (Cathedral Cave, MO) & HRV
%     % (Harvard) -- please note that the input order is not important:
%     newpf=extractpair(pf,'CCM','HRV');
%
%     % Pull profiles for station CMB of the Berkeley network:
%     newpf=extractpair(pf,{'BK' 'CMB'});
%
%    See also: SLOWDECAYPAIRS, PLOT_CMB_MEASUREMENTS, PLOT_CMB_PDF

%     Version History:
%        Feb.  2, 2011 - initial version
%        Mar. 30, 2011 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 30, 2011 at 13:35 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% check profile struct
reqfields={'gcdist','azwidth','slow','slowerr','decay','decayerr',...
    'cslow','cslowerr','cdecay','cdecayerr','cluster','kname','st','ev',...
    'delaz','synthetics','earthmodel','corrections','corrcoef','freq',...
    'phase','runname','dirname','time'};
if(~isstruct(pf) || any(~isfield(pf,reqfields)))
    error('seizmo:extractpair:badInput',...
        ['PF must be a struct with the fields:\n' ...
        sprintf('''%s'' ',reqfields{:}) '!']);
end

% check stn info
if(isstring(stn1)); stn1=cellstr(stn1); end
if(~iscellstr(stn1) || numel(stn1)>4)
    error('seizmo:extractpair:badInput',...
        ['STN1 must be one of the following:\n' ...
        '''NAME'',{''NET'' ''NAME''},{''NET'' ''NAME'' ''STREAM''}']);
end
if(nargin==3)
    if(isstring(stn2)); stn2=cellstr(stn2); end
    if(~iscellstr(stn2) || numel(stn2)>4)
        error('seizmo:extractpair:badInput',...
            ['STN2 must be one of the following:\n' ...
            '''NAME'',{''NET'' ''NAME''},{''NET'' ''NAME'' ''STREAM''}']);
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
    warning('seizmo:extractpair:no2StationFound',...
        'No 2-station profiles given!');
    return;
end

% combine kname info (each pair is a new page)
kname=upper(reshape([kname{:}],2,4,[]));
npf=size(kname,3);

% find matches for first station
match11=false(npf,1); match12=false(npf,1);
stn1=upper(stn1); n1=numel(stn1);
switch n1
    case 1 % name
        match11(:)=strcmpi(stn1,kname(1,2,:));
        match12(:)=strcmpi(stn1,kname(2,2,:));
    case 2 % net name
        for i=1:npf
            match11(i)=isequal(stn1,kname(1,1:2,i));
            match12(i)=isequal(stn1,kname(2,1:2,i));
        end
    case 3 % net name hole
        for i=1:npf
            match11(i)=isequal(stn1,kname(1,1:3,i));
            match12(i)=isequal(stn1,kname(2,1:3,i));
        end
    case 4 % net name hole cmp
        for i=1:npf
            match11(i)=isequal(stn1,kname(1,1:4,i));
            match12(i)=isequal(stn1,kname(2,1:4,i));
        end
end

% find matches for station 2
if(nargin==3)
    match21=false(npf,1); match22=false(npf,1);
    stn2=upper(stn2); n2=numel(stn2);
    switch n2
        case 1 % name
            match21(:)=strcmpi(stn2,kname(1,2,:));
            match22(:)=strcmpi(stn2,kname(2,2,:));
        case 2 % net name
            for i=1:npf
                match21(i)=isequal(stn2,kname(1,1:2,i));
                match22(i)=isequal(stn2,kname(2,1:2,i));
            end
        case 3 % net name hole
            for i=1:npf
                match21(i)=isequal(stn2,kname(1,1:3,i));
                match22(i)=isequal(stn2,kname(2,1:3,i));
            end
        case 4 % net name hole cmp
            for i=1:npf
                match21(i)=isequal(stn2,kname(1,1:4,i));
                match22(i)=isequal(stn2,kname(2,1:4,i));
            end
    end
end

if(nargin==3)
    % find pair matches
    pf=pf((match11 & match22) | (match12 & match21));
else
    % single station matches
    pf=pf(match11 | match12);
end

end
