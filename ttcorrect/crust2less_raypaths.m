function [paths]=crust2less_raypaths(paths)
%CRUST2LESS_RAYPATHS    Removes ray path above Crust2.0 moho
%
%    Usage:    paths=crust2less_raypaths(paths)
%
%    Description: PATHS=CRUST2LESS_RAYPATHS(PATHS) will remove the
%     portions of the raypaths that are above the Crust2.0 moho.  To avoid
%     plotting functions connecting across segments where the crust was
%     removed, NaNs are inserted.  Segments that cross the moho are
%     adjusted so that the last point is at the moho boundary.
%
%    Notes:
%     - Currently this isn't very smart.  It will have trouble if it
%       encounters a crust2.0 sidewall.  A warning is issued and no
%       adjustment is done (actually the crossing segment is removed).
%
%    Examples:
%     This is mainly meant for MANCOR.  Once plotting is up, I'll add an
%     example here.
%
%    See also: GETRAYPATHS, GET_UPSWING_RAYPATHS, MANCOR, TAUPPATH

%     Version History:
%        May  31, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  31, 2010 at 20:05 GMT

% todo:
% - test for lat/lon less paths
%
% We don't handle the case where the path intersects the side of a
% crustal wall.  How to do it?  If moho is below both points then
% we have this case.  Then we need to find intersection with wall (not
% really easy).  Currently we just drop the crossing segment.
%
%  ___   
%   x |    ___
%     | x  |
%     |____|

% check nargin
error(nargchk(1,1,nargin));

% input should be a tauppath struct
test=tauppath('ph','P','deg',10);
if(~isstruct(paths) || any(~ismember(fieldnames(paths),fieldnames(test))))
    error('seizmo:crust2less_raypaths:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
elseif(any(~ismember(fieldnames(paths(1).path),fieldnames(test(1).path))))
    error('seizmo:crust2less_raypaths:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
end

% loop over each path
for i=1:numel(paths)
    % check for lat/lon-less paths
    if(isequal(unique(paths(i).path.latitude),0) ...
            && isequal(unique(paths(i).path.longitude),0))
        warning('seizmo:crust2less_raypaths:noLatLon',...
            'PATHS appears to be missing lat/lon info!');
    end
    
    % find points below moho
    moho=getc2moho(paths(i).path.latitude,paths(i).path.longitude);
    below=moho<paths(i).path.depth;
    npts=numel(moho);
    
    % now get start/end indices of sections
    trans=diff([false; below; false]);
    s=find(trans==1);
    e=find(trans==-1)-1;
    
    % handle nothing/everything below moho
    if(isempty(s))
        % set all to empty
        paths(i).path.rayparameter=[];
        paths(i).path.time=[];
        paths(i).path.distance=[];
        paths(i).path.depth=[];
        paths(i).path.latitude=[];
        paths(i).path.longitude=[];
        continue;
    elseif(isequal(s,1) && isequal(e,npts))
        % everything below so go to next
        continue;
    end
    
    % extend sections to include segments crossing the moho
    s=s-1; s(s<1)=1;
    e=e+1; e(e>npts)=npts;
    ns=numel(s);
    
    % loop over sections below
    newrp=cell(ns,1); newt=newrp; newdist=newrp;
    newdep=newrp; newlat=newrp; newlon=newrp;
    for j=1:ns
        % extract section
        newrp{j}=paths(i).path.rayparameter(s(j):e(j));
        newt{j}=paths(i).path.time(s(j):e(j));
        newdist{j}=paths(i).path.distance(s(j):e(j));
        newdep{j}=paths(i).path.depth(s(j):e(j));
        newlat{j}=paths(i).path.latitude(s(j):e(j));
        newlon{j}=paths(i).path.longitude(s(j):e(j));
        
        % fix segments crossing moho
        if(~below(s(j)))
            % check if both are above the moho (sidewall condition)
            % - warning b/c I don't handle this properly yet
            % - this is not the definitive test for sidewall but is
            %   the case that can break things
            if(all(newdep{j}(1:2)<=moho(s(j))))
                warning('seizmo:crust2less_raypaths:sidewall',...
                    'Hit a crust side wall: path %d, section %d start!',...
                    i,j);
                
                % just chop off first point and move on
                newdep{j}(1)=[];
                newrp{j}(1)=[];
                newt{j}(1)=[];
                newdist{j}(1)=[];
                newlat{j}(1)=[];
                newlon{j}(1)=[];
            else
                % get fraction of crossing segment below moho
                % - assumes moho is between two points
                frac=abs((moho(s(j))-newdep{j}(2))...
                    /(diff(newdep{j}(1:2))));

                % shift values of point above moho to the moho
                % - assumes values behave linearly
                % - also assumes we don't hit a sidewall of the crust
                % - lat/lon are handled specially to avoid issues
                newdep{j}(1)=moho(s(j));
                newrp{j}(1)=newrp{j}(2)-frac*diff(newrp{j}(1:2));
                newt{j}(1)=newt{j}(2)-frac*diff(newt{j}(1:2));
                newdist{j}(1)=newdist{j}(2)-frac*diff(newdist{j}(1:2));
                [dist,az]=sphericalinv(...
                    paths(i).path.latitude(2),paths(i).path.longitude(2),...
                    paths(i).path.latitude(1),paths(i).path.longitude(1));
                [paths(i).path.latitude(1),paths(i).path.longitude(1)]=...
                    sphericalfwd(paths(i).path.latitude(2),...
                    paths(i).path.longitude(2),dist*frac,az);
            end
        end
        if(~below(e(j)))
            % check if both are above the moho (sidewall condition)
            % - warning b/c I don't handle this properly yet
            % - this is not the definitive test for sidewall but is
            %   the case that can break things
            if(all(newdep{j}(end-1:end)<=moho(e(j))))
                warning('seizmo:crust2less_raypaths:sidewall',...
                    'Hit a crust side wall: path %d, section %d end!',...
                    i,j);
                
                % just chop off last point and move on
                newdep{j}(end)=[];
                newrp{j}(end)=[];
                newt{j}(end)=[];
                newdist{j}(end)=[];
                newlat{j}(end)=[];
                newlon{j}(end)=[];
            else
                % get fraction of crossing segment below moho
                % - assumes moho is between two points
                frac=abs((moho(e(j))-newdep{j}(end-1))...
                    /(diff(newdep{j}(end-1:end))));
                
                % shift values of point above moho to the moho
                % - assumes values behave linearly
                % - also assumes we don't hit a sidewall of the crust
                % - lat/lon are handled specially to avoid issues
                newdep{j}(end)=moho(e(j));
                newrp{j}(end)=newrp{j}(end-1)...
                    +frac*diff(newrp{j}(end-1:end));
                newt{j}(end)=newt{j}(end-1)...
                    +frac*diff(newt{j}(end-1:end));
                newdist{j}(end)=newdist{j}(end-1)...
                    +frac*diff(newdist{j}(end-1:end));
                [dist,az]=sphericalinv(...
                    paths(i).path.latitude(end-1),...
                    paths(i).path.longitude(end-1),...
                    paths(i).path.latitude(end),...
                    paths(i).path.longitude(end));
                [paths(i).path.latitude(end),...
                    paths(i).path.longitude(end)]=...
                    sphericalfwd(paths(i).path.latitude(end-1),...
                    paths(i).path.longitude(end-1),dist*frac,az);
            end
        end
        
        % add nan to end of each section (except for the last one)
        if(j~=ns)
            newrp{j}=[newrp{j}; nan];
            newt{j}=[newt{j}; nan];
            newdist{j}=[newdist{j}; nan];
            newdep{j}=[newdep{j}; nan];
            newlat{j}=[newlat{j}; nan];
            newlon{j}=[newlon{j}; nan];
        end
    end
    
    % and concatenate new sections together (with nans between)
    paths(i).path.rayparameter=cat(1,newrp{:});
    paths(i).path.time=cat(1,newt{:});
    paths(i).path.distance=cat(1,newdist{:});
    paths(i).path.depth=cat(1,newdep{:});
    paths(i).path.latitude=cat(1,newlat{:});
    paths(i).path.longitude=cat(1,newlon{:});
end

end
