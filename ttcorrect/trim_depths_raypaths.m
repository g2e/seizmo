function [paths]=trim_depths_raypaths(paths,range)
%TRIM_DEPTHS_RAYPATHS    Removes sections of raypaths outside depth range
%
%    Usage:    paths=trim_depths_raypaths(paths,range)
%
%    Description:
%     PATHS=TRIM_DEPTHS_RAYPATHS(PATHS,RANGE) will remove the portions of
%     the raypaths in struct PATHS that are outside the depth range RANGE.
%     PATHS should be a struct matching the layout as output from TAUPPATH
%     (GETRAYPATHS conforms).  PATHS is adjusted so that points outside the
%     depth range are removed and segments of the path that cross the depth
%     boundary are adjusted so that they are on the boundary using linear
%     interpolation.  Individual path sections within the range are
%     separated by NaNs as to avoid inappropriate connection during
%     plotting (or in MANCOR).
%
%    Notes:
%
%    Examples:
%     % Only show raypath sections within the core:
%     paths=tauppath('ev',[5 129],'st',[41 -1]);
%     paths=trim_depths_raypaths(paths,[2891 6371]);
%     plotraypaths(paths);
%
%    See also: CRUSTLESS_RAYPATHS, GETRAYPATHS, PLOTRAYPATHS, TAUPPATH,
%              INSERT_DEPTHS_IN_RAYPATHS, EXTRACT_UPSWING_RAYPATHS, MANCOR

%     Version History:
%        June  2, 2010 - initial version
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
    error('seizmo:trim_depths_raypaths:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
elseif(any(~ismember(fieldnames(paths(1).path),fieldnames(test(1).path))))
    error('seizmo:trim_depths_raypaths:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
elseif(~isreal(range) || numel(range)~=2 || ~issorted(range))
    error('seizmo:trim_depths_raypaths:badRANGE',...
        'RANGE must be a real-valued vector [top bottom] in km!');
end

% specific locations or not?
% - this fails if some are and some are not
if(any(~ismember({'latitude' 'longitude'},fieldnames(paths(1).path))))
    specific=false;
else % not specific
    specific=true;
end

% loop over each path
for i=1:numel(paths)
    % find points inside range
    in=paths(i).path.depth>=range(1) & paths(i).path.depth<=range(2);
    npts=numel(paths(i).path.depth);
    
    % now get start/end indices of sections
    % - note that sections that hop over the range are not included!
    trans=diff([false; in; false]);
    s=find(trans==1);
    e=find(trans==-1)-1;
    
    % handle nothing/everything in range
    if(isempty(s))
        % set all to empty
        paths(i).path.time=[];
        paths(i).path.distance=[];
        paths(i).path.depth=[];
        if(specific)
            paths(i).path.latitude=[];
            paths(i).path.longitude=[];
        end
        continue;
    elseif(isequal(s,1) && isequal(e,npts))
        % everything below so go to next
        continue;
    end
    
    % extend sections to include segments crossing boundaries
    s=s-1; s(s<1)=1;
    e=e+1; e(e>npts)=npts;
    ns=numel(s);
    
    % loop over sections
    newt=cell(ns,1); newdist=newt;
    newdep=newt; newlat=newt; newlon=newt;
    for j=1:ns
        % extract section
        newt{j}=paths(i).path.time(s(j):e(j));
        newdist{j}=paths(i).path.distance(s(j):e(j));
        newdep{j}=paths(i).path.depth(s(j):e(j));
        if(specific)
            newlat{j}=paths(i).path.latitude(s(j):e(j));
            newlon{j}=paths(i).path.longitude(s(j):e(j));
        end
        
        % fix boundary crossing segments
        if(~in(s(j)) && ~isnan(newdep{j}(1)))
            % which boundary are we crossing?
            if(newdep{j}(1)<range(1))
                % top
                frac=abs((range(1)-newdep{j}(2))...
                    /(diff(newdep{j}(1:2))));
                newdep{j}(1)=range(1);
            else
                % bottom
                frac=abs((range(2)-newdep{j}(2))...
                    /(diff(newdep{j}(1:2))));
                newdep{j}(1)=range(2);
            end

            % shift values of point to the boundary
            % - assumes values behave linearly
            % - lat/lon are handled specially to avoid issues
            newt{j}(1)=newt{j}(2)-frac*diff(newt{j}(1:2));
            newdist{j}(1)=newdist{j}(2)-frac*diff(newdist{j}(1:2));
            if(specific)
                [dist,az]=sphericalinv(paths(i).path.latitude(2),...
                    paths(i).path.longitude(2),...
                    paths(i).path.latitude(1),paths(i).path.longitude(1));
                [paths(i).path.latitude(1),paths(i).path.longitude(1)]=...
                    sphericalfwd(paths(i).path.latitude(2),...
                    paths(i).path.longitude(2),dist*frac,az);
            end
        elseif(isnan(newdep{j}(1)))
            % just chop off first point and move on
            newdep{j}(1)=[];
            newt{j}(1)=[];
            newdist{j}(1)=[];
            if(specific)
                newlat{j}(1)=[];
                newlon{j}(1)=[];
            end
        end
        if(~in(e(j)) && ~isnan(newdep{j}(end)))
            % which boundary are we crossing?
            if(newdep{j}(end)<range(1))
                % top
                frac=abs((range(1)-newdep{j}(end-1))...
                    /(diff(newdep{j}(end-1:end))));
                newdep{j}(end)=range(1);
            else
                % bottom
                frac=abs((range(2)-newdep{j}(end-1))...
                    /(diff(newdep{j}(end-1:end))));
                newdep{j}(end)=range(2);
            end

            % shift values of point to the boundary
            % - assumes values behave linearly
            % - lat/lon are handled specially to avoid issues
            newt{j}(end)=newt{j}(end-1)+frac*diff(newt{j}(end-1:end));
            newdist{j}(end)=newdist{j}(end-1)...
                +frac*diff(newdist{j}(end-1:end));
            if(specific)
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
        elseif(isnan(newdep{j}(end)))
            % just chop off first point and move on
            newdep{j}(end)=[];
            newt{j}(end)=[];
            newdist{j}(end)=[];
            if(specific)
                newlat{j}(end)=[];
                newlon{j}(end)=[];
            end
        end
        
        % add nan to end of each section (except for the last one)
        if(j~=ns)
            newt{j}=[newt{j}; nan];
            newdist{j}=[newdist{j}; nan];
            newdep{j}=[newdep{j}; nan];
            if(specific)
                newlat{j}=[newlat{j}; nan];
                newlon{j}=[newlon{j}; nan];
            end
        end
    end
    
    % and concatenate new sections together (with nans between)
    paths(i).path.time=cat(1,newt{:});
    paths(i).path.distance=cat(1,newdist{:});
    paths(i).path.depth=cat(1,newdep{:});
    if(specific)
        paths(i).path.latitude=cat(1,newlat{:});
        paths(i).path.longitude=cat(1,newlon{:});
    end
end

end
