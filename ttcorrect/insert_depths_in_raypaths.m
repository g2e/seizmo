function [paths]=insert_depths_in_raypaths(paths,depths)
%INSERT_DEPTHS_IN_RAYPATHS    Add specific depths to raypath points
%
%    Usage:    paths=insert_depths_in_raypaths(paths,depths)
%
%    Description: PATHS=INSERT_DEPTHS_IN_RAYPATHS(PATHS,DEPTHS) adds points
%     to the raypaths in struct PATHS at the depths indicated in DEPTHS.
%     PATHS should conform to the output from TAUPPATH.  DEPTHS is expected
%     to be in units of km.  Notes that depths are only added if they are
%     crossed by the raypaths.
%
%    Notes:
%
%    Examples:
%     Insert some points at all PREM discontinuities:
%      mod=prem:
%      depths=mod.depth(find(~diff(mod.depth)))';
%      paths=insert_depths_in_raypaths(paths,depths);
%
%     Effectively upsampling your raypath:
%      paths=insert_depths_in_raypaths(paths,0:5:2890);
%
%    See also: GETRAYPATHS, TRIM_DEPTHS_RAYPATHS, TAUPPATH,
%              CRUST2LESS_RAYPATHS, GET_UPSWING_RAYPATHS

%     Version History:
%        June  3, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  3, 2010 at 23:15 GMT

% todo:

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
                paths(i).path.rayparameter=...
                    [paths(i).path.rayparameter(1:cnt+j);
                     paths(i).path.rayparameter(cnt+j)+frac*(diff(paths(i).path.rayparameter(cnt+j:cnt+j+1)));
                     paths(i).path.rayparameter(cnt+j+1:end)];
                paths(i).path.time=...
                    [paths(i).path.time(1:cnt+j);
                     paths(i).path.time(cnt+j)+frac*(diff(paths(i).path.time(cnt+j:cnt+j+1)));
                     paths(i).path.time(cnt+j+1:end)];
                paths(i).path.distance=...
                    [paths(i).path.distance(1:cnt+j);
                     paths(i).path.distance(cnt+j)+frac*(diff(paths(i).path.distance(cnt+j:cnt+j+1)));
                     paths(i).path.distance(cnt+j+1:end)];
                [dist,az]=sphericalinv(...
                    paths(i).path.latitude(cnt+j),paths(i).path.longitude(cnt+j),...
                    paths(i).path.latitude(cnt+j+1),paths(i).path.longitude(cnt+j+1));
                [lat,lon]=sphericalfwd(...
                    paths(i).path.latitude(cnt+j),paths(i).path.longitude(cnt+j),dist*frac,az);
                paths(i).path.latitude=...
                    [paths(i).path.latitude(1:cnt+j);
                     lat;
                     paths(i).path.latitude(cnt+j+1:end)];
                paths(i).path.longitude=...
                    [paths(i).path.longitude(1:cnt+j);
                     lon;
                     paths(i).path.longitude(cnt+j+1:end)];
                
                % increment counter for subsequent tests
                cnt=cnt+1;
            end
        end
    end
end
