function [paths]=crustless_raypaths(paths,model)
%CRUSTLESS_RAYPATHS    Removes ray path above moho
%
%    Usage:    paths=crustless_raypaths(paths)
%              paths=crustless_raypaths(paths,model)
%
%    Description:
%     PATHS=CRUSTLESS_RAYPATHS(PATHS) will remove the portions of the
%     raypaths that are above the moho.  To avoid plotting functions
%     connecting across segments where the crust was removed, NaNs are
%     inserted to isolate the mantle raypath sections from each other.
%     Segments that cross the moho are adjusted so that the segments end at
%     the moho boundary.  Linear interpolation is used to estimate these
%     positions.
%
%     PATHS=CRUSTLESS_RAYPATHS(PATHS,MODEL) selects the crustal model.  The
%     default model and available models are defined by GETCRUST.
%
%    Notes:
%     - This is not a very smart function.  It will have trouble if it
%       encounters a crustal sidewall.  One way this condition occurs is
%       when a raypath segment crosses a block boundary AND both ends of
%       the segment are above the moho depth of the block that the segment
%       ends in but below the moho that the segment starts in.  To work
%       around this case the crossing segment is removed.  A warning is
%       issued if seizmodebug=true so you know if this happens.
%
%    Examples:
%     % Plot some paths without the crustal portions:
%     paths=tauppath('ev',[5 129],'st',[41 -1]);
%     plotraypaths(crustless_raypaths(paths));
%
%    See also: GETRAYPATHS, EXTRACT_UPSWING_RAYPATHS, MANCOR, TAUPPATH,
%              TRIM_DEPTHS_RAYPATHS, INSERT_DEPTHS_IN_RAYPATHS,
%              PLOTRAYPATHS

%     Version History:
%        May  31, 2010 - initial version
%        June  3, 2010 - updated example, handle nans
%        Aug.  8, 2010 - doc update
%        Feb. 27, 2012 - update for tauppath changes
%        Aug.  6, 2012 - warn only if seizmodebug=true
%        Jan. 23, 2014 - update for crustal model changes, model option
%        Jan. 31, 2014 - fixed error when no model specified
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 31, 2014 at 02:45 GMT

% todo:
% We don't handle the case where the path intersects the side of a
% crustal wall.  How to do it?  If one of the moho depths is below both
% points then we have this case.  Then we need to find intersection with
% wall (not really easy and this is not even correct to match crust and
% mantle corrections).  Currently we just drop the crossing segment.
%
% Both of these segments would be dropped:
%  ___   
%   x |    ___
%     | x  |
%     |____| x

% check nargin
error(nargchk(1,2,nargin));

% input should be a tauppath struct
test=tauppath('ph','P','ev',[0 0],'st',[0 10]);
if(~isstruct(paths) || any(~ismember(fieldnames(paths),fieldnames(test))))
    error('seizmo:crustless_raypaths:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
elseif(any(~ismember(fieldnames(paths(1).path),fieldnames(test(1).path))))
    error('seizmo:crustless_raypaths:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
elseif(any(~ismember({'latitude' 'longitude'},fieldnames(paths(1).path))))
    error('seizmo:crustless_raypaths:badStruct',...
        'Latitude & Longitude are required path fields!');
end

% default/check model
if(nargin<2); model=[]; end
if(~isempty(model) ...
        && (~ischar(model) || ndims(model)~=2 || size(model,1)~=1))
    error('seizmo:crustless_raypaths:badInput',...
        'MODEL must be a string!');
end

% hide warnings unless debugging is on
debug=seizmodebug;

% loop over each path
for i=1:numel(paths)
    % check for lat/lon-less paths
    if(isequal(unique(paths(i).path.latitude),0) ...
            && isequal(unique(paths(i).path.longitude),0))
        warning('seizmo:crustless_raypaths:noLatLon',...
            'PATHS appears to be missing lat/lon info!');
    end
    
    % check for nans
    nn=~isnan(paths(i).path.depth);
    npts=numel(paths(i).path.depth);
    
    % find points below moho
    mod=getcrust(paths(i).path.latitude(nn),...
        paths(i).path.longitude(nn),model);
    moho=nan(npts,1);
    moho(nn)=-mod.top(:,9);
    below=moho<paths(i).path.depth;
    
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
    newt=cell(ns,1); newdist=newt;
    newdep=newt; newlat=newt; newlon=newt;
    for j=1:ns
        % extract section
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
            if(isnan(newdep{j}(1)))
                % just chop off first point and move on
                newdep{j}(1)=[];
                newt{j}(1)=[];
                newdist{j}(1)=[];
                newlat{j}(1)=[];
                newlon{j}(1)=[];
            elseif(all(newdep{j}(1:2)<=moho(s(j))))
                if(debug)
                    warning('seizmo:crustless_raypaths:sidewall',...
                        ['Hit a crust side wall: ' ...
                        'path %d, section %d start!'],i,j);
                end
                
                % just chop off first point and move on
                newdep{j}(1)=[];
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
            if(isnan(newdep{j}(end)))
                % just chop off last point and move on
                newdep{j}(end)=[];
                newt{j}(end)=[];
                newdist{j}(end)=[];
                newlat{j}(end)=[];
                newlon{j}(end)=[];
            elseif(all(newdep{j}(end-1:end)<=moho(e(j))))
                if(debug)
                    warning('seizmo:crustless_raypaths:sidewall',...
                        ['Hit a crust side wall: ' ...
                        'path %d, section %d end!'],i,j);
                end
                
                % just chop off last point and move on
                newdep{j}(end)=[];
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
            newt{j}=[newt{j}; nan];
            newdist{j}=[newdist{j}; nan];
            newdep{j}=[newdep{j}; nan];
            newlat{j}=[newlat{j}; nan];
            newlon{j}=[newlon{j}; nan];
        end
    end
    
    % and concatenate new sections together (with nans between)
    paths(i).path.time=cat(1,newt{:});
    paths(i).path.distance=cat(1,newdist{:});
    paths(i).path.depth=cat(1,newdep{:});
    paths(i).path.latitude=cat(1,newlat{:});
    paths(i).path.longitude=cat(1,newlon{:});
end

end
