function [db,idx]=db_update(db,varargin)
%DB_UPDATE    Update db and/or return index for field with the value

try
    % get values already in db
    values=getsubfield(db,varargin{1:end-1});
    nv=numel(values); idx=[];

    % if values is a cell string array we need to unwrap that
    if(iscellstr(values))
        % look for a match
        for i=1:nv
            if(isequalwithequalnans(varargin{end},values{i}))
                idx=i;
            end
        end
    else % not a cell string so do not unwrap cells
        for i=1:nv
            if(isequalwithequalnans(varargin{end},values(i)))
                idx=i;
            end
        end
    end

    % new value, so update by adding to the end
    if(isempty(idx))
        idx=nv+1;
        db=setfield(db,varargin{1:end-1},idx,varargin{end});
    end
catch
    idx=1;
    db=setfield(db,varargin{1:end-1},idx,varargin{end});
end

end
