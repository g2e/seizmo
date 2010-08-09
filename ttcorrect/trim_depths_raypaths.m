function [paths]=trim_depths_raypaths(paths,range)
%TRIM_DEPTHS_RAYPATHS    Removes sections of raypaths outside depth range
%
%    Usage:    paths=trim_depths_raypaths(paths,range)
%
%    Description: PATHS=TRIM_DEPTHS_RAYPATHS(PATHS,RANGE) will remove the
%     portions of the raypaths in struct PATHS that are outside the depth
%     range RANGE.  PATHS should be a struct matching the layout as output
%     from TAUPPATH (GETRAYPATH conforms).  PATHS is adjusted so that
%     points outside the depth range are removed and segments of the path
%     that cross the depth boundary are adjusted so that they are on the
%     boundary using linear interpolation.  Individual path sections within
%     the range are separated by NaNs as to avoid inappropriate connection
%     during plotting (or in MANCOR).
%
%    Notes:
%
%    Examples:
%     Only show raypath sections within the core:
%      paths=tauppath('ev',[5 129],'st',[41 -1]);
%      paths=trim_depths_raypaths(paths,[2891 6371]);
%      plotraypaths(paths);
%
%    See also: CRUST2LESS_RAYPATHS, GETRAYPATHS, PLOTRAYPATHS,
%              INSERT_DEPTHS_IN_RAYPATHS, EXTRACT_UPSWING_RAYPATHS, MANCOR

%     Version History:
%        June  2, 2010 - initial version
%        Aug.  8, 2010 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  8, 2010 at 23:15 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% input should be a tauppath struct
test=tauppath('ph','P','deg',10);
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
    
    % extend sections to include segments crossing boundaries
    s=s-1; s(s<1)=1;
    e=e+1; e(e>npts)=npts;
    ns=numel(s);
    
    % loop over sections
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
            newrp{j}(1)=newrp{j}(2)-frac*diff(newrp{j}(1:2));
            newt{j}(1)=newt{j}(2)-frac*diff(newt{j}(1:2));
            newdist{j}(1)=newdist{j}(2)-frac*diff(newdist{j}(1:2));
            [dist,az]=sphericalinv(...
                paths(i).path.latitude(2),paths(i).path.longitude(2),...
                paths(i).path.latitude(1),paths(i).path.longitude(1));
            [paths(i).path.latitude(1),paths(i).path.longitude(1)]=...
                sphericalfwd(paths(i).path.latitude(2),...
                paths(i).path.longitude(2),dist*frac,az);
        elseif(isnan(newdep{j}(1)))
            % just chop off first point and move on
            newdep{j}(1)=[];
            newrp{j}(1)=[];
            newt{j}(1)=[];
            newdist{j}(1)=[];
            newlat{j}(1)=[];
            newlon{j}(1)=[];
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
            newrp{j}(end)=newrp{j}(end-1)+frac*diff(newrp{j}(end-1:end));
            newt{j}(end)=newt{j}(end-1)+frac*diff(newt{j}(end-1:end));
            newdist{j}(end)=newdist{j}(end-1)+frac*diff(newdist{j}(end-1:end));
            [dist,az]=sphericalinv(...
                paths(i).path.latitude(end-1),paths(i).path.longitude(end-1),...
                paths(i).path.latitude(end),paths(i).path.longitude(end));
            [paths(i).path.latitude(end),paths(i).path.longitude(end)]=...
                sphericalfwd(paths(i).path.latitude(end-1),...
                paths(i).path.longitude(end-1),dist*frac,az);
        elseif(isnan(newdep{j}(end)))
            % just chop off first point and move on
            newdep{j}(end)=[];
            newrp{j}(end)=[];
            newt{j}(end)=[];
            newdist{j}(end)=[];
            newlat{j}(end)=[];
            newlon{j}(end)=[];
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
