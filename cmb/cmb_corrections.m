function [varargout]=cmb_corrections(varargin)
%CMB_CORRECTIONS    Gets travel time and amplitude corrections
%
%    Usage:    results=cmb_corrections(results)
%              corrections=cmb_corrections(phase,data)
%
%    Description:
%     RESULTS=CMB_CORRECTIONS(RESULTS) addes travel time and amplitude
%     correction info to the struct RESULTS (created by CMB_1ST_PASS,
%     CMB_OUTLIERS, or CMB_2ND_PASS).  Please note that this is called by
%     CMB_1ST_PASS already, so there is no need to recompute the
%     corrections unless you have changed the corrections in some manner.
%     The corrections are stored in the RESULTS struct under the field
%     ".corrections".
%
%     CORRECTIONS=CMB_CORRECTIONS(PHASE,DATA) returns a struct CORRECTIONS
%     containing travel time and amplitude corrections for the seismic
%     phase given by PHASE contained in the seismic records in DATA.  PHASE
%     is a string such as 'Pdiff' or 'Sdiff'.  DATA is a SEIZMO data struct
%     (use "help seizmo" for more info).
%
%    Notes:
%     - Currently only supports phases:
%        Pdiff  &  Sdiff
%
%    Examples:
%     % Get corrections for Pdiff for your dataset:
%     corrections=cmb_corrections('Pdiff',data);
%
%    See also: CMB_1ST_PASS, CMB_OUTLIERS, CMB_2ND_PASS, SLOWDECAYPAIRS,
%              SLOWDECAYPROFILES, MAP_CMB_PROFILES, PREP_CMB_DATA,
%              TEST_SCALED_MANTLE_CORRECTIONS

%     Version History:
%        Oct.  7, 2010 - initial version
%        Oct. 11, 2010 - fix bug in populating geom spreading corrections
%        Dec. 29, 2010 - added docs, 2nd usage type simplifies my life
%        Jan. 12, 2011 - fixed several bugs from new usage format
%        Jan. 29, 2011 - use check_cmb_results
%        Mar. 10, 2011 - skip all mantle corrections but HMSL, save paths,
%                        radiation pattern amplitudes
%        Jan. 23, 2014 - minor fix for a rename of called function
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 23, 2014 at 15:25 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check results struct if single input
% - then extract appropriate info and set output flag
if(nargin==1)
    % checking
    error(check_cmb_results(varargin{1}));
    
    % recall using appropriate fields
    varargout{1}=varargin{1};
    for i=1:numel(varargin{1})
        if(~isempty(varargin{1}(i).useralign))
            varargout{1}(i).corrections=...
                cmb_corrections(varargin{1}(i).phase,...
                varargin{1}(i).useralign.data);
        else
            varargout{1}(i).corrections=[];
        end
    end
    return;
end

% check phase
valid={'Pdiff' 'SHdiff' 'SVdiff'};
if(~isstring(varargin{1}) || ~ismember(varargin{1},valid))
    error('seizmo:cmb_corrections:badPhase',...
        ['PHASE must be one of the following:\n' ...
        sprintf('''%s'' ',valid{:}) '!']);
end

% check data
error(seizmocheck(varargin{2}));

% necessary header info
% ev - evla evlo evel evdp
% st - stla stlo stel stdp
% delaz - gcarc az baz dist
[ev,delaz,st,kevnm]=getheader(varargin{2},'ev','delaz','st','kevnm');

% get cmt from kevnm
kevnm=unique(kevnm);
if(numel(kevnm)~=1 || strcmpi(kevnm,'NaN'))
    error('seizmo:cmb_corrections:badData',...
        'Cannot find CMT info for this data!');
end
cmt=findcmt('name',char(kevnm));

% convert meters to kilometers
ev(:,3:4)=ev(:,3:4)/1000;
st(:,3:4)=st(:,3:4)/1000;

% operation depends on phase
switch varargin{1}
    case 'Pdiff'
        % get ellipticity corrections
        corrections.ellcor=ellcor(ev(:,1),ev(:,4),delaz(:,1),delaz(:,2),'Pdiff');
        
        % get crustal corrections
        rayp=4.42802574759071;
        corrections.crucor.prem=crucor(st(:,1),st(:,2),rayp,'P',...
            'elev',st(:,3),'hole',st(:,4),'refmod','prem');
        %corrections.crucor.ak135=crucor(st(:,1),st(:,2),rayp,'P',...
        %    'elev',st(:,3),'hole',st(:,4),'refmod','ak135');
        %corrections.crucor.iasp91=crucor(st(:,1),st(:,2),rayp,'P',...
        %    'elev',st(:,3),'hole',st(:,4),'refmod','iasp91');
        
        % get Sdiff raypaths
        corrections.paths=getraypaths('P,Pdiff','prem',ev(:,1),ev(:,2),ev(:,4),st(:,1),st(:,2));
        
        % remove crust from paths
        corrections.paths=crustless_raypaths(corrections.paths);
        
        % upswing paths
        % - using 500km above CMB as the cutoff
        corrections.uppaths=extract_upswing_raypaths(corrections.paths,2890-500);
        
        % mantle corrections
        corrections.mancor.hmsl06p.full=mancor(corrections.paths,'hmsl06p');
        corrections.mancor.hmsl06p.upswing=mancor(corrections.uppaths,'hmsl06p');
        %corrections.mancor.mitp08.full=mancor(paths,'mit-p08');
        %corrections.mancor.mitp08.upswing=mancor(uppaths,'mit-p08');
        %corrections.mancor.dz04.full=mancor(paths,'dz04');
        %corrections.mancor.dz04.upswing=mancor(uppaths,'dz04');
        %corrections.mancor.prip05.full=mancor(paths,'pri-p05');
        %corrections.mancor.prip05.upswing=mancor(uppaths,'pri-p05');
    case {'SHdiff' 'SVdiff'}
        % get ellipticity corrections
        corrections.ellcor=ellcor(ev(:,1),ev(:,4),delaz(:,1),delaz(:,2),'Sdiff');
        
        % get crustal corrections
        rayp=8.36067454903639;
        corrections.crucor.prem=crucor(st(:,1),st(:,2),rayp,'S',...
            'elev',st(:,3),'hole',st(:,4),'refmod','prem');
        %corrections.crucor.ak135=crucor(st(:,1),st(:,2),rayp,'S',...
        %    'elev',st(:,3),'hole',st(:,4),'refmod','ak135');
        %corrections.crucor.iasp91=crucor(st(:,1),st(:,2),rayp,'S',...
        %    'elev',st(:,3),'hole',st(:,4),'refmod','iasp91');
        
        % get Sdiff raypaths
        corrections.paths=getraypaths('S,Sdiff','prem',ev(:,1),ev(:,2),ev(:,4),st(:,1),st(:,2));
        
        % remove crust from paths
        corrections.paths=crustless_raypaths(corrections.paths);
        
        % upswing paths
        % - using 500km above CMB as the cutoff
        corrections.uppaths=extract_upswing_raypaths(corrections.paths,2890-500);
        
        % mantle corrections
        corrections.mancor.hmsl06s.upswing=mancor(corrections.uppaths,'hmsl06s');
        corrections.mancor.hmsl06s.full=mancor(corrections.paths,'hmsl06s');
        %corrections.mancor.s20rts.upswing=mancor(uppaths,'s20rts');
        %corrections.mancor.s20rts.full=mancor(paths,'s20rts');
        %corrections.mancor.saw24b16.upswing=mancor(uppaths,'saw24b16');
        %corrections.mancor.saw24b16.full=mancor(paths,'saw24b16');
        %corrections.mancor.sb4l18.upswing=mancor(uppaths,'sb4l18');
        %corrections.mancor.sb4l18.full=mancor(paths,'sb4l18');
        %corrections.mancor.tx2007.upswing=mancor(uppaths,'tx2007');
        %corrections.mancor.tx2007.full=mancor(paths,'tx2007');
        %corrections.mancor.pris05.full=mancor(paths,'pri-s05');
        %corrections.mancor.pris05.upswing=mancor(uppaths,'pri-s05');
end

% get geometrical spreading corrections
corrections.geomsprcor=geomsprcor(delaz(:,1));

% radiation pattern amplitude corrections
model=prem('depth',unique(ev(:,4)));
inc=rayp2inc([corrections.paths.rayparameter]',...
    model.(['v' lower(varargin{1}(1))]),6371-unique(ev(:,4)));
corrections.radpatcor=radpat(cmt,inc,delaz(:,2),varargin{1}(1:end-4));

% output
varargout{1}=corrections;

end
