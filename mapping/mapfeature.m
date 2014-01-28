function []=mapfeature(varargin)
%MAPFEATURE    Map geographic feature
%
%    Usage:    mapfeature('feature')
%              mapfeature({'feature1' 'feature2' ...})
%              mapfeature(features,opts)
%              mapfeature(features,opts,'prop1',val1,'prop2',val2,...)
%              mapfeature(ax,...)
%
%    Description:
%     MAPFEATURE('FEATURE') maps the feature given by string 'FEATURE'.
%     Currently the feature set is:
%      casz                 - central african shear zone
%      cob                  - continent-ocean crustal boundary
%      congo_craton         - Congo craton
%      cvl                  - Cameroon volcanic line
%      fitton_benue         - Fitton's Benue trough
%      fitton_cvl           - Fitton's Cameroon volcanic line
%      fraczone             - fracture zones
%      hotspots             - hotspots
%      impacts              - confirmed impact craters
%      isochrons            - oceanic crust isochrons
%      josplateau           - jos plateau granites
%      lips                 - large igneous provinces (incomplete)
%      maglines             - magnetic lineaments
%      pb2002_boundaries    - Peter Bird plate boundaries
%      pb2002_orogens       - Peter Bird orogenic regions
%      pb2002_plates        - Peter Bird plates
%      ridge                - ridge plate boundaries (not pb2002)
%      seamounts            - seamounts
%      sleep_hotspots       - Norm Sleep's hotspot locations
%      transform            - transform plate boundaries (not pb2002)
%      trench               - convergent plate boundaries (not pb2002)
%      volcanoes            - volcanoes
%      wars                 - west african rift system
%      xridge               - extinct ridges
%     Most features have shortcut names - see the code to find them.
%
%     MAPFEATURE({'FEATURE1' 'FEATURE2' ...}) maps multiple features.  The
%     features must be strings in a cell array.  They are plotted in the
%     order given.
%
%     MAPFEATURE(FEATURES,OPTS) specifies how polygonal features are
%     mapped.  The choices are 'lines', 'patch', & 'hatch'.  Use 'lines' to
%     map the polygon outlines, 'patch' to draw them as patches, and
%     'hatch' to draw outlines & a hatched internal area.  You may also use
%     a cell array triplet {'style' angle step} to implicitly specify a
%     hatch and adjust its hatching.  See M_HATCH for details.  Use a cell
%     array if you wish to specify a different style for each feature in
%     FEATURES.
%
%     MAPFEATURE(FEATURES,OPTS,'PROP1',VAL1,'PROP2',VAL2,...) passes
%     property/value pairs to M_LINE, M_PATCH, or M_HATCH.  Which in turn
%     pass the properties on to LINE & PATCH.
%
%     MAPFEATURE(AX,...) draws on the M_Map map given by axes handle AX.
%     Please note that M_Map cannot handle drawing amongst multiple map
%     instances sequentially (eg you must draw 1 map at a time).
%
%    Notes:
%     - All features are tagged based on the dataset name.  See the source
%       code below if you need help finding the tag for a feature.
%     - LIPS is only a partial subset.  If you want to help clean it up see
%       the TODO section in the code below.
%     - COB is incomplete and needs help.  If you want to help, use
%       the crustal thickness of a crustal model as a guide.
%     - Dataset sources: see FEATURE_REFERENCES in See Also section
%
%    Examples:
%     % Map plate boundaries, seamounts & hotspots:
%     mapfeature({'boundaries' 'seamounts' 'hotspots'})
%
%     % Map plates as patches and diffuse plate boundaries in cross-hatch:
%     mapfeature({'plates' 'orogens'},{'p' {'cross' 45 1}})
%
%     % Make a map focused on central africa and add some features:
%     ax=mmap('po',{'lat',[-20 20],'lon',[-20 20]},'go',{'box','fancy'});
%     mapfeature(ax,{'boundaries' 'isochrons' 'fz' 'seamounts' 'lips' ...
%                    'wars' 'cvl' 'jos' 'cc' 'impacts' 'volcanoes' ...
%                    'hotspots'});
%
%    See also: MMAP, FEATURE_REFERENCES, RAISEFANCY

%     Version History:
%        Feb.  9, 2011 - initial version
%        Feb. 16, 2012 - fix bug in parsing hatch triplets
%        Jan. 23, 2014 - minor doc update
%        Jan. 27, 2014 - use axparse instead of axescheck for octave
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 27, 2014 at 23:00 GMT

% todo:
% - fix the polygon insanity of some features
%   - all polygons should proceed clockwise (use sph_poly_area)
%   - lips is just crazy
%     - current "clean" set is all closed, simple polygons oriented cw
%     - final set rules:
%       - remove points & lines (ie less than 3 unique points = 248)
%       - remove line features (what are they?)
%         - this is tough as some are parts of a broken polygon
%       - fix simple self-intersecting (this is the most common)
%         - 2 crossing segments separated by 1 segment
%       - fix those with overlapping segments (at least 1)
%       - make sure all proceed clockwise
%       - remove duplicate points (unclose polygons)
% - europe cob is botched
% - can we split peter bird boundaries into ridge, trench, transform?
%   - ridge/transform are not differentiated
%   - convergent/sub-left/sub-right are

% extract axes handle
[ax,varargin]=axparse(varargin{:});
nargs=numel(varargin);

% check nargin
error(nargchk(1,inf,nargs));
if(nargs>2 && mod(nargs,2))
    error('seizmo:mapfeature:badInput',...
        'All Properties must be paired with a value!');
end

% extract feature & option
if(nargs>0); feature=varargin{1}; varargin(1)=[]; end
if(nargs>1); opt=varargin{1}; varargin(1)=[]; end

% require feature to be a string or cellstr
if(isstring(feature)); feature=cellstr(feature); end
if(~iscellstr(feature))
    error('seizmo:mapfeature:badInput',...
        'FEATURE must be a string or cell array of strings!');
end
nf=numel(feature);

% default opt (1 empty array)
if(nargs==1 || isempty(opt)); opt={[]}; end

% check opt
hopt=cell(nf,1);
if(isstring(opt)); opt=cellstr(opt); end
if(iscell(opt))
    % encapsulate {'...' ## ##}
    if(numel(opt)==3 && isstring(opt{1}) && isnumeric(opt{2}) ...
            && isnumeric(opt{3}))
        opt={opt};
    end
    
    % require scalar or nf elements
    if(~any(numel(opt)==[1 nf]))
        error('seizmo:mapfeature:badInput',...
            'OPT must be a single entry or as many entries as features!');
    end
    
    % check each element
    for i=1:numel(opt)
        % require option to be either:
        % [], 'lines', 'patches', 'hatches' or {'...' ## ##}
        if(isempty(opt{i})); continue; end
        if(isstring(opt{i}))
            if(strmatch(lower(opt{i}),'lines')); opt{i}='lines';
            elseif(strmatch(lower(opt{i}),'patches')); opt{i}='patch';
            elseif(strmatch(lower(opt{i}),'hatches')); opt{i}='hatch';
            else
                error('seizmo:mapfeature:badInput',...
                    ['OPT string must be either ''LINES'' ''PATCHES'' ' ...
                    'or ''HATCHES''!']);
            end
            continue;
        end
        if(numel(opt{i})==3 && isstring(opt{i}{1}) ...
                && isnumeric(opt{i}{2}) && isnumeric(opt{i}{3}))
            hopt{i}=opt{i};
            opt{i}='hatch';
            continue;
        end
        error('seizmo:mapfeature:badInput',...
            'OPT specified incorrectly!');
    end
    
    % expand scalar opt
    if(isscalar(opt)); opt=opt(ones(nf,1)); end
else
    error('seizmo:mapfeature:badInput',...
        'OPT specified incorrectly!');
end

% check axis (new figure if none)
if(isempty(ax))
    % new map (full globe)
    ax=mmap;
elseif(isscalar(ax) && isreal(ax) ...
        && ishandle(ax) && strcmp('axes',get(ax,'type')))
    axes(ax);
else
    error('seizmo:mapfeature:badInput',...
        'AX must be a valid axes handle!');
end

% require map
global MAP_PROJECTION MAP_VAR_LIST
if(isempty(MAP_PROJECTION))
    error('seizmo:mapfeature:badInput',...
        'AX must have been created using M_Map!');
end

% check/set hold
held=ishold(ax);
hold(ax,'on');

% loop over features
for a=1:numel(feature)
    switch lower(feature{a})
        case {'central africa shear zone' 'casz'}
            % black lines of width 2
            f=load('casz');
            for b=1:numel(f.casz)
                f.casz(b).longitude=...
                    unwrap(f.casz(b).longitude*pi/180)*180/pi;
                m_line(f.casz(b).longitude,f.casz(b).latitude,...
                    'color','k','linewidth',2,'tag','casz','parent',ax,...
                    varargin{:});
            end
        case {'continent ocean boundary' 'cob'}
            % black lines of width 1
            f=load('cob');
            for b=1:numel(f.cob)
                f.cob(b).longitude=...
                    unwrap(f.cob(b).longitude*pi/180)*180/pi;
                m_line(f.cob(b).longitude,f.cob(b).latitude,...
                    'color','k','linewidth',1,'tag','cob','parent',ax,...
                    varargin{:});
            end
        case {'congo craton' 'congo_craton' 'cc'}
            % patch, line, or hatched polygons
            f=load('congo_craton');
            if(isempty(opt{a}))
                opt{a}='hatch';
                hopt{a}=[];
            end
            for b=1:numel(f.congo_craton)
                f.congo_craton(b).longitude=...
                    unwrap(f.congo_craton(b).longitude*pi/180)*180/pi;
                switch opt{a}
                    case 'lines'
                        % blue dashed line of width 2
                        m_line(f.congo_craton(b).longitude,...
                            f.congo_craton(b).latitude,...
                            'color','b','linewidth',2,'linestyle','--',...
                            'tag','congo_craton','parent',ax,varargin{:});
                    case 'patch'
                        % blue patch with black dashed outline of width 2
                        m_patch(f.congo_craton(b).longitude,...
                            f.congo_craton(b).latitude,...
                            'b','linewidth',2,'linestyle','--',...
                            'tag','congo_craton','parent',ax,varargin{:});
                    case 'hatch'
                        % blue dashed line of width 2 surrounding hatch of
                        % width 1 oriented at 45deg spaced at 1 pixel
                        if(isempty(hopt{a}))
                            hopt{a}={'single' 45 1};
                        end
                        m_hatch(f.congo_craton(b).longitude,...
                            f.congo_craton(b).latitude,hopt{a}{:},...
                            'color','b','linewidth',1,'linestyle','-',...
                            'tag','congo_craton','parent',ax,varargin{:});
                        m_line(f.congo_craton(b).longitude,...
                            f.congo_craton(b).latitude,...
                            'color',name2rgb('b')/2,'linewidth',2,...
                            'linestyle','--','tag','congo_craton',...
                            'parent',ax,varargin{:});
                end
            end
        case {'cameroon volcanic line' 'cvl'}
            % patches, lines, or hatched polygons
            f=load('cvl');
            if(isempty(opt{a}))
                opt{a}='patch';
                hopt{a}=[];
            end
            for b=1:numel(f.cvl)
                f.cvl(b).longitude=...
                    unwrap(f.cvl(b).longitude*pi/180)*180/pi;
                switch opt{a}
                    case 'lines'
                        m_line(f.cvl(b).longitude,f.cvl(b).latitude,...
                            'color','y','linewidth',1,...
                            'tag','cvl','parent',ax,varargin{:});
                    case 'patch'
                        m_patch(f.cvl(b).longitude,f.cvl(b).latitude,...
                            'y','linewidth',1,...
                            'tag','cvl','parent',ax,varargin{:});
                    case 'hatch'
                        if(isempty(hopt{a}))
                            hopt{a}={'cross' 45 1};
                        end
                        m_hatch(f.cvl(b).longitude,f.cvl(b).latitude,...
                            hopt{a}{:},'color','y','linewidth',1,...
                            'tag','cvl','parent',ax,varargin{:});
                        m_line(f.cvl(b).longitude,f.cvl(b).latitude,...
                            'color',name2rgb('y')/3,'linewidth',2,...
                            'tag','cvl','parent',ax,varargin{:});
                end
            end
        case {'fitton benue trough' 'fitton_benue' 'benue'}
            % patch, line, or hatched polygons
            f=load('fitton_benue');
            if(isempty(opt{a}))
                opt{a}='patch';
                hopt{a}=[];
            end
            for b=1:numel(f.fitton_benue)
                f.fitton_benue(b).longitude=...
                    unwrap(f.fitton_benue(b).longitude*pi/180)*180/pi;
                switch opt{a}
                    case 'lines'
                        m_line(f.fitton_benue(b).longitude,...
                            f.fitton_benue(b).latitude,...
                            'color',name2rgb('orange'),'linewidth',1,...
                            'tag','fitton_benue','parent',ax,varargin{:});
                    case 'patch'
                        m_patch(f.fitton_benue(b).longitude,...
                            f.fitton_benue(b).latitude,...
                            name2rgb('orange'),'linewidth',1,...
                            'tag','fitton_benue','parent',ax,varargin{:});
                    case 'hatch'
                        if(isempty(hopt{a}))
                            hopt{a}={'cross' 15 1};
                        end
                        m_hatch(f.fitton_benue(b).longitude,...
                            f.fitton_benue(b).latitude,hopt{a}{:},...
                            'color',name2rgb('orange'),'linewidth',1,...
                            'tag','fitton_benue','parent',ax,varargin{:});
                        m_line(f.fitton_benue(b).longitude,...
                            f.fitton_benue(b).latitude,...
                            'color',name2rgb('orange')/2,'linewidth',2,...
                            'tag','fitton_benue','parent',ax,varargin{:});
                end
            end
        case {'fitton cameroon volcanic line' 'fitton cvl' 'fitton_cvl'}
            % patches, lines, or hatched polygons
            f=load('fitton_cvl');
            if(isempty(opt{a}))
                opt{a}='patch';
                hopt{a}=[];
            end
            for b=1:numel(f.fitton_cvl)
                f.fitton_cvl(b).longitude=...
                    unwrap(f.fitton_cvl(b).longitude*pi/180)*180/pi;
                switch opt{a}
                    case 'lines'
                        m_line(f.fitton_cvl(b).longitude,...
                            f.fitton_cvl(b).latitude,'color','y',...
                            'linewidth',1,'tag','fitton_cvl',...
                            'parent',ax,varargin{:});
                    case 'patch'
                        m_patch(f.fitton_cvl(b).longitude,...
                            f.fitton_cvl(b).latitude,'y','linewidth',1,...
                            'tag','fitton_cvl','parent',ax,varargin{:});
                    case 'hatch'
                        if(isempty(hopt{a}))
                            hopt{a}={'cross' 45 1};
                        end
                        m_hatch(f.fitton_cvl(b).longitude,...
                            f.fitton_cvl(b).latitude,hopt{a}{:},...
                            'color','y','linewidth',1,...
                            'tag','fitton_cvl','parent',ax,varargin{:});
                        m_line(f.fitton_cvl(b).longitude,...
                            f.fitton_cvl(b).latitude,...
                            'color',name2rgb('y')/3,...
                            'linewidth',2,'tag','fitton_cvl',...
                            'parent',ax,varargin{:});
                end
            end
        case {'fracture zones' 'fractures' 'fraczones' 'fraczone' 'fz'}
            % black lines
            f=load('fraczone');
            for b=1:numel(f.fraczone)
                f.fraczone(b).longitude=...
                    unwrap(f.fraczone(b).longitude*pi/180)*180/pi;
                m_line(f.fraczone(b).longitude,f.fraczone(b).latitude,...
                    'color','k','linewidth',1,'tag','fraczone',...
                    'parent',ax,varargin{:});
            end
        case 'hotspots'
            % orange & red circles
            f=load('hotspots');
            lonrange=diff(MAP_VAR_LIST.longs);
            pixel=30/(lonrange/180);
            if(pixel>100); pixel=100; end
            lons=[f.hotspots.longitude].';
            lats=[f.hotspots.latitude].';
            while(any(abs(lons-mean(MAP_VAR_LIST.longs))>180))
                lons(lons<mean(MAP_VAR_LIST.longs)-180)=...
                    lons(lons<mean(MAP_VAR_LIST.longs)-180)+360;
                lons(lons>mean(MAP_VAR_LIST.longs)+180)=...
                    lons(lons>mean(MAP_VAR_LIST.longs)+180)-360;
            end
            m_scatter(lons,lats,pixel,name2rgb('orange'),'o',...
                'markerfacecolor','r','parent',ax,...
                'tag','hotspots',varargin{:});
        case 'impacts'
            % violet asterisks
            f=load('impacts');
            pixel=[f.impacts.diameter].';
            lons=[f.impacts.longitude].';
            lats=[f.impacts.latitude].';
            while(any(abs(lons-mean(MAP_VAR_LIST.longs))>180))
                lons(lons<mean(MAP_VAR_LIST.longs)-180)=...
                    lons(lons<mean(MAP_VAR_LIST.longs)-180)+360;
                lons(lons>mean(MAP_VAR_LIST.longs)+180)=...
                    lons(lons>mean(MAP_VAR_LIST.longs)+180)-360;
            end
            m_scatter(lons,lats,pixel,name2rgb('violet'),'*',...
                'parent',ax,'tag','impacts',varargin{:});
        case 'isochrons'
            % hsv-colored lines (blue is old, red is young)
            f=load('isochrons');
            niso=numel(f.isochrons);
            ages=[f.isochrons.age]';
            cmap=[interp1q([0; 1],[5/6; 0],1-ages/200) ones(niso,2)];
            cmap(cmap<0)=0;
            cmap(cmap>1)=1;
            cmap=hsv2rgb(cmap);
            for b=1:niso
                f.isochrons(b).longitude=...
                    unwrap(f.isochrons(b).longitude*pi/180)*180/pi;
                m_line(f.isochrons(b).longitude,f.isochrons(b).latitude,...
                    'color',cmap(b,:),'linewidth',1,'tag','isochrons',...
                    'parent',ax,varargin{:});
            end
        case {'jos plateau granites' 'josplateau' 'jos'}
            % patches, lines, or hatched polygons
            f=load('josplateau');
            if(isempty(opt{a}))
                opt{a}='patch';
                hopt{a}=[];
            end
            for b=1:numel(f.josplateau)
                f.josplateau(b).longitude=...
                    unwrap(f.josplateau(b).longitude*pi/180)*180/pi;
                switch opt{a}
                    case 'lines'
                        m_line(f.josplateau(b).longitude,...
                            f.josplateau(b).latitude,...
                            'color','r','linewidth',1,...
                            'tag','josplateau','parent',ax,varargin{:});
                    case 'patch'
                        m_patch(f.josplateau(b).longitude,...
                            f.josplateau(b).latitude,...
                            'r','linewidth',1,...
                            'tag','josplateau','parent',ax,varargin{:});
                    case 'hatch'
                        if(isempty(hopt{a}))
                            hopt{a}={'cross' 45 1};
                        end
                        m_hatch(f.josplateau(b).longitude,...
                            f.josplateau(b).latitude,hopt{a}{:},...
                            'color','r','linewidth',1,...
                            'tag','josplateau','parent',ax,varargin{:});
                        m_line(f.josplateau(b).longitude,...
                            f.josplateau(b).latitude,...
                            'color',[.5 0 0],'linewidth',2,...
                            'tag','josplateau','parent',ax,varargin{:});
                end
            end
        case {'large igneous provinces' 'lips'}
            % patches, lines, or hatched polygons
            f=load('lips_clean');
            if(isempty(opt{a}))
                opt{a}='patch';
                hopt{a}=[];
            end
            for b=1:numel(f.lips)
                f.lips(b).longitude=...
                    unwrap(f.lips(b).longitude*pi/180)*180/pi;
                switch opt{a}
                    case 'lines'
                        m_line(f.lips(b).longitude,f.lips(b).latitude,...
                            'color','y','linewidth',1,...
                            'tag','lips','parent',ax,varargin{:});
                    case 'patch'
                        m_patch(f.lips(b).longitude,f.lips(b).latitude,...
                            'y','linewidth',1,...
                            'tag','lips','parent',ax,varargin{:});
                    case 'hatch'
                        if(isempty(hopt{a}))
                            hopt{a}={'cross' 45 1};
                        end
                        m_hatch(f.lips(b).longitude,f.lips(b).latitude,...
                            hopt{a}{:},'color','y','linewidth',1,...
                            'tag','lips','parent',ax,varargin{:});
                        m_line(f.lips(b).longitude,f.lips(b).latitude,...
                            'color',name2rgb('y')/3,'linewidth',2,...
                            'tag','lips','parent',ax,varargin{:});
                end
            end
        case {'magnetic lineaments' 'magnetic lines' 'maglines'}
            % lines
            f=load('maglines');
            for b=1:numel(f.maglines)
                f.maglines(b).longitude=...
                    unwrap(f.maglines(b).longitude*pi/180)*180/pi;
                m_line(f.maglines(b).longitude,f.maglines(b).latitude,...
                    'color','g','linewidth',1,'tag','maglines',...
                    'parent',ax,varargin{:});
            end
        case {'pb2002_boundaries' 'peter bird plate boundaries' ...
                'plate boundaries' 'boundaries'}
            % lines
            f=load('pb2002_boundaries');
            for b=1:numel(f.pb2002_boundaries)
                f.pb2002_boundaries(b).longitude=...
                    unwrap(f.pb2002_boundaries(b).longitude*pi/180)*180/pi;
                m_line(f.pb2002_boundaries(b).longitude,...
                    f.pb2002_boundaries(b).latitude,...
                    'color','r','linewidth',1,'tag','pb2002_boundaries',...
                    'parent',ax,varargin{:});
            end
        case {'pb2002_orogens' 'peter bird orogens' 'orogens' ...
                'orogenic zones' 'diffuse plate boundaries'}
            % patches, lines, hatched polygons
            f=load('pb2002_orogens');
            if(isempty(opt{a}))
                opt{a}='hatch';
                hopt{a}=[];
            end
            for b=1:numel(f.pb2002_orogens)
                f.pb2002_orogens(b).longitude=...
                    unwrap(f.pb2002_orogens(b).longitude*pi/180)*180/pi;
                switch opt{a}
                    case 'lines'
                        % solid red outlines
                        m_line(f.pb2002_orogens(b).longitude,...
                            f.pb2002_orogens(b).latitude,...
                            'color','r','linewidth',1,...
                            'tag','pb2002_orogens',...
                            'parent',ax,varargin{:});
                    case 'patch'
                        % black outlined red patches
                        m_patch(f.pb2002_orogens(b).longitude,...
                            f.pb2002_orogens(b).latitude,...
                            'r','linewidth',1,...
                            'tag','pb2002_orogens',...
                            'parent',ax,varargin{:});
                    case 'hatch'
                        % red hatched areas
                        if(isempty(hopt{a}))
                            hopt{a}={'single' 45 1};
                        end
                        m_hatch(f.pb2002_orogens(b).longitude,...
                            f.pb2002_orogens(b).latitude,hopt{a}{:},...
                            'color','r','linewidth',1,...
                            'tag','pb2002_orogens',...
                            'parent',ax,varargin{:});
                        m_line(f.pb2002_orogens(b).longitude,...
                            f.pb2002_orogens(b).latitude,...
                            'color',name2rgb('r')/2,'linewidth',2,...
                            'tag','pb2002_orogens',...
                            'parent',ax,varargin{:});
                end
            end
        case {'pb2002_plates' 'peter bird plates' 'plates'}
            % patches, lines, hatched polygons
            f=load('pb2002_plates');
            if(isempty(opt{a}))
                opt{a}='patch';
                hopt{a}=[];
            end
            np=numel(f.pb2002_plates);
            cmap=hsv(np);
            rand('seed',666);
            cmap=cmap(randperm(np),:);
            rand('seed',now);
            for b=1:numel(f.pb2002_plates)
                f.pb2002_plates(b).longitude=...
                    unwrap(f.pb2002_plates(b).longitude*pi/180)*180/pi;
                switch opt{a}
                    case 'lines'
                        % solid red outlines
                        m_line(f.pb2002_plates(b).longitude,...
                            f.pb2002_plates(b).latitude,...
                            'color','r','linewidth',1,...
                            'tag','pb2002_plates','parent',ax,varargin{:});
                    case 'patch'
                        % hsv colored patches (black outlines)
                        m_patch(f.pb2002_plates(b).longitude,...
                            f.pb2002_plates(b).latitude,...
                            cmap(b,:),'linewidth',1,...
                            'tag','pb2002_plates','parent',ax,varargin{:});
                    case 'hatch'
                        % hsv-colored hatched areas
                        if(isempty(hopt{a}))
                            hopt{a}={'cross' 45 1};
                        end
                        m_hatch(f.pb2002_plates(b).longitude,...
                            f.pb2002_plates(b).latitude,hopt{a}{:},...
                            'color',cmap(b,:),'linewidth',1,...
                            'tag','pb2002_plates','parent',ax,varargin{:});
                        m_line(f.pb2002_plates(b).longitude,...
                            f.pb2002_plates(b).latitude,...
                            'color','r','linewidth',2,...
                            'tag','pb2002_plates','parent',ax,varargin{:});
                end
            end
        case {'simple plate boundaries' 'simple boundaries'}
            % solid red lines
            f=load('plates');
            for b=1:numel(f.plates)
                f.plates(b).longitude=...
                    unwrap(f.plates(b).longitude*pi/180)*180/pi;
                m_line(f.plates(b).longitude,f.plates(b).latitude,...
                    'color','r','linewidth',1,...
                    'tag','plates','parent',ax,varargin{:});
            end
        case {'ridge' 'ridges'}
            % solid red lines
            f=load('ridge');
            for b=1:numel(f.ridge)
                f.ridge(b).longitude=...
                    unwrap(f.ridge(b).longitude*pi/180)*180/pi;
                m_line(f.ridge(b).longitude,f.ridge(b).latitude,...
                    'color','r','linewidth',1,...
                    'tag','ridge','parent',ax,varargin{:});
            end
        case 'seamounts'
            % aquamarine dots
            f=load('seamounts');
            lonrange=diff(MAP_VAR_LIST.longs);
            pixel=10/(lonrange/360);
            if(pixel>100); pixel=100; end
            for b=1:numel(f.seamounts)
                lons=[f.seamounts(b).longitude].';
                lats=[f.seamounts(b).latitude].';
                while(any(abs(lons-mean(MAP_VAR_LIST.longs))>180))
                    lons(lons<mean(MAP_VAR_LIST.longs)-180)=...
                        lons(lons<mean(MAP_VAR_LIST.longs)-180)+360;
                    lons(lons>mean(MAP_VAR_LIST.longs)+180)=...
                        lons(lons>mean(MAP_VAR_LIST.longs)+180)-360;
                end
                m_scatter(lons,lats,pixel,name2rgb('aquamarine'),'.',...
                    'parent',ax,'tag','seamounts',varargin{:});
            end
        case {'sleep hotspots' 'sleep_hotspots'}
            % red-filled purple circles
            f=load('sleep_hotspots');
            lonrange=diff(MAP_VAR_LIST.longs);
            pixel=30/(lonrange/180);
            if(pixel>100); pixel=100; end
            lons=[f.sleep_hotspots.longitude].';
            lats=[f.sleep_hotspots.latitude].';
            while(any(abs(lons-mean(MAP_VAR_LIST.longs))>180))
                lons(lons<mean(MAP_VAR_LIST.longs)-180)=...
                    lons(lons<mean(MAP_VAR_LIST.longs)-180)+360;
                lons(lons>mean(MAP_VAR_LIST.longs)+180)=...
                    lons(lons>mean(MAP_VAR_LIST.longs)+180)-360;
            end
            m_scatter(lons,lats,pixel,name2rgb('violet'),'o',...
                'markerfacecolor','r','parent',ax,...
                'tag','sleep_hotspots',varargin{:});
        case {'transform' 'transforms' 'transform faults'}
            % solid green lines
            f=load('transform');
            for b=1:numel(f.transform)
                f.transform(b).longitude=...
                    unwrap(f.transform(b).longitude*pi/180)*180/pi;
                m_line(f.transform(b).longitude,f.transform(b).latitude,...
                    'color','g','linewidth',1,...
                    'tag','transform','parent',ax,varargin{:});
            end
        case {'trench' 'trenches' 'subduction' 'subduction zones'}
            % solid blue lines
            f=load('trench');
            for b=1:numel(f.trench)
                f.trench(b).longitude=...
                    unwrap(f.trench(b).longitude*pi/180)*180/pi;
                m_line(f.trench(b).longitude,f.trench(b).latitude,...
                    'color','b','linewidth',1,...
                    'tag','trench','parent',ax,varargin{:});
            end
        case {'volcanoes' 'volcanos'}
            % brown-filled upward-pointing red triangles
            f=load('volcanoes');
            lonrange=diff(MAP_VAR_LIST.longs);
            pixel=10/(lonrange/360);
            if(pixel>100); pixel=100; end
            lons=[f.volcanoes.longitude].';
            lats=[f.volcanoes.latitude].';
            while(any(abs(lons-mean(MAP_VAR_LIST.longs))>180))
                lons(lons<mean(MAP_VAR_LIST.longs)-180)=...
                    lons(lons<mean(MAP_VAR_LIST.longs)-180)+360;
                lons(lons>mean(MAP_VAR_LIST.longs)+180)=...
                    lons(lons>mean(MAP_VAR_LIST.longs)+180)-360;
            end
            m_scatter(lons,lats,pixel,[139 69 19]/255,'^',...
                'markerfacecolor','r','parent',ax,...
                'tag','volcanoes',varargin{:});
        case {'west african rift system' 'wars'}
            % patches, lines, or hatched polygons
            f=load('wars');
            if(isempty(opt{a}))
                opt{a}='patch';
                hopt{a}=[];
            end
            for b=1:numel(f.wars)
                f.wars(b).longitude=...
                    unwrap(f.wars(b).longitude*pi/180)*180/pi;
                switch opt{a}
                    case 'lines'
                        % orange solid outline
                        m_line(f.wars(b).longitude,f.wars(b).latitude,...
                            'color',name2rgb('orange'),'linewidth',1,...
                            'tag','wars','parent',ax,varargin{:});
                    case 'patch'
                        % black outlined orange patches
                        m_patch(f.wars(b).longitude,f.wars(b).latitude,...
                            name2rgb('orange'),'linewidth',1,...
                            'tag','wars','parent',ax,varargin{:});
                    case 'hatch'
                        % orange hatched areas
                        if(isempty(hopt{a}))
                            hopt{a}={'cross' 45 1};
                        end
                        m_hatch(f.wars(b).longitude,f.wars(b).latitude,...
                            hopt{a}{:},'color',name2rgb('orange'),...
                            'linewidth',1,'tag','wars',...
                            'parent',ax,varargin{:});
                        m_line(f.wars(b).longitude,f.wars(b).latitude,...
                            'color',name2rgb('orange')/2,'linewidth',2,...
                            'tag','wars','parent',ax,varargin{:});
                end
            end
        case {'extinct ridges' 'xridges' 'xridge'}
            % red dashed lines over gray solid lines
            f=load('xridge');
            for b=1:numel(f.xridge)
                f.xridge(b).longitude=...
                    unwrap(f.xridge(b).longitude*pi/180)*180/pi;
                m_line(f.xridge(b).longitude,f.xridge(b).latitude,...
                    'color',[0.5 0.5 0.5],'linewidth',2,...
                    'tag','xridge','parent',ax,varargin{:});
                m_line(f.xridge(b).longitude,f.xridge(b).latitude,...
                    'color','r','linewidth',1,'linestyle','-.',...
                    'tag','xridge','parent',ax,varargin{:});
            end
        otherwise
            warning('seizmo:mapfeature:badInput',...
                'Cannot map unknown feature: %s',feature{i});
    end
end

% set hold back
if(~held); hold(ax,'off'); end

end
