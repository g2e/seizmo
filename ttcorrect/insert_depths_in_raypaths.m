function [paths]=insert_depths_in_raypaths(paths,depths)
%INSERT_DEPTHS_IN_RAYPATHS    Add specific depths to raypath points
%
%    Usage:    paths=insert_depths_in_raypaths(paths,depths)
%
%    Description:
%     PATHS=INSERT_DEPTHS_IN_RAYPATHS(PATHS,DEPTHS) adds points to the
%     raypaths in struct PATHS at the depths indicated in DEPTHS using
%     linear interpolation.  PATHS should conform to the output from
%     TAUPPATH.  DEPTHS is expected to be in units of km.  Notes that
%     depths are only added if they are crossed by the raypaths.
%
%    Notes:
%     - SPHERICALFWD & SPHERICALINV are used to make sure the raypaths do
%       not have issues at the poles.
%
%    Examples:
%     % Insert some points at all PREM discontinuities:
%     mod=prem:
%     depths=mod.depth(find(~diff(mod.depth)))';
%     paths=insert_depths_in_raypaths(paths,depths);
%
%     % Effectively upsampling your raypath:
%     paths=insert_depths_in_raypaths(paths,0:5:2890);
%
%    See also: GETRAYPATHS, TRIM_DEPTHS_RAYPATHS, TAUPPATH, PLOTRAYPATHS,
%              CRUSTLESS_RAYPATHS, EXTRACT_UPSWING_RAYPATHS, MANCOR

%     Version History:
%        June  3, 2010 - initial version
%        Aug.  8, 2010 - doc update
%        Feb. 27, 2012 - update for tauppath changes
%        Jan. 23, 2014 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 23, 2014 at 02:45 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% input should be a tauppath struct
test=tauppath('ph','P','ev',[0 0],'st',[0 10]);
if(~isstruct(paths) || any(~ismember(fieldnames(paths),fieldnames(test))))
    error('seizmo:insert_depths_in_raypaths:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
elseif(any(~ismember(fieldnames(paths(1).path),fieldnames(test(1).path))))
    error('seizmo:insert_depths_in_raypaths:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
end

% specific locations or not?
% - this fails if some are and some are not
if(any(~ismember({'latitude' 'longitude'},fieldnames(paths(1).path))))
    specific=false;
else % not specific
    specific=true;
end

% number of raypaths
nrp=numel(paths);
ndep=numel(depths);

% sort depths
depths=sort(depths(:));

% loop over raypaths
for i=1:nrp
    % number of segments
    nseg=numel(paths(i).path.depth)-1;
    
    % loop over each segment
    cnt=0;
    for j=1:nseg
        % incrementing over depths depends on direction of segment
        dd=paths(i).path.depth(cnt+j+1)-paths(i).path.depth(cnt+j);
        if(dd<0)
            % upgoing so go deep to shallow
            idx=ndep:-1:1;
        elseif(dd>0)
            % downgoing so go shallow to deep
            idx=1:ndep;
        else
            % same depth so no need to look
            continue;
        end
        
        % loop over depths
        for k=idx
            if(dd>0)
                mind=paths(i).path.depth(cnt+j);
                maxd=paths(i).path.depth(cnt+j+1);
            else
                mind=paths(i).path.depth(cnt+j+1);
                maxd=paths(i).path.depth(cnt+j);
            end
            
            % is depth in segment?
            if(depths(k)<maxd && depths(k)>mind)
                % YES
                
                % get fraction of depth difference of depth from first
                frac=abs((depths(k)-paths(i).path.depth(cnt+j))/dd);
                
                % insert depth
                paths(i).path.depth=...
                    [paths(i).path.depth(1:cnt+j);
                     depths(k);
                     paths(i).path.depth(cnt+j+1:end)];
                paths(i).path.time=...
                    [paths(i).path.time(1:cnt+j);
                     paths(i).path.time(cnt+j)...
                        +frac*(diff(paths(i).path.time(cnt+j:cnt+j+1)));
                     paths(i).path.time(cnt+j+1:end)];
                paths(i).path.distance=...
                    [paths(i).path.distance(1:cnt+j);
                     paths(i).path.distance(cnt+j)...
                        +frac*(diff(paths(i).path.distance(cnt+j:cnt+j+1)));
                     paths(i).path.distance(cnt+j+1:end)];
                if(specific)
                    [dist,az]=sphericalinv(...
                        paths(i).path.latitude(cnt+j),...
                        paths(i).path.longitude(cnt+j),...
                        paths(i).path.latitude(cnt+j+1),...
                        paths(i).path.longitude(cnt+j+1));
                    [lat,lon]=sphericalfwd(...
                        paths(i).path.latitude(cnt+j),...
                        paths(i).path.longitude(cnt+j),dist*frac,az);
                    paths(i).path.latitude=...
                        [paths(i).path.latitude(1:cnt+j);
                         lat;
                         paths(i).path.latitude(cnt+j+1:end)];
                    paths(i).path.longitude=...
                        [paths(i).path.longitude(1:cnt+j);
                         lon;
                         paths(i).path.longitude(cnt+j+1:end)];
                end
                
                % increment counter for subsequent tests
                cnt=cnt+1;
            end
        end
    end
end
