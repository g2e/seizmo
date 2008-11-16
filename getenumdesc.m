function [varargout]=getenumdesc(data,varargin)
%GETENUMDESC    Get enum string description from SEIZMO header enum field
%
%    Description: GETENUMDESC(DATA,FIELD) returns a cellstring array filled 
%     with the description strings associated with the enum field FIELD
%     stored in the SEIZMO structure DATA.  This is more useful/readible
%     than the magic number (enumerator integer) returned with GETHEADER.
%
%     GETENUMDESC(DATA,FIELD1,...,FIELDN) returns a cellstring array for
%     each field supplied.
%
%    Notes:
%     - Numeric header fields not defined to be enum fields can be used as 
%       if they were enum fields with GETENUMDESC.  This gives the user the 
%       ability to have more enum fields if needed.  Character fields are
%       NOT able to be treated as enum fields.
%     - Nonexistent header fields will return 'Unknown Enum Field'.
%     - Enum values that aren't defined in the enumerator list are tagged
%       as 'Unknown Enum Value' unless they are the defined UNDEFINED
%       value in which case they are tagged as 'Undefined Enum Field'. This
%       avoids conflict with enum fields called 'unknown' and 'undefined'.
%
%    Tested on: Matlab r2007b
%
%    Usage: cellstr=getenumdesc(data,'field')
%           [cellstr1,...,cellstrN]=getenumdesc(data,'field1',...,'fieldN')
%
%    Examples:
%     To check if all records are timeseries data:
%      if(all(strcmp(getenumdesc(data,'iftype'),'Time Series File')))
%          disp('timeseries dataset')
%      end
%
%     Interpret resp0 as an enum field:
%      my_enum_desc=getenumdesc(data,'resp0')   
%
%    See also: getheader, getlgc, getenumid

%     Version History:
%        Feb. 23, 2008 - initial version
%        Feb. 28, 2008 - minor code clean
%        Feb. 29, 2008 - handle undefined
%        Mar.  4, 2008 - minor doc update
%        June 13, 2008 - doc update, compat fixes, handle values out of
%                        range or non-whole, avoid name conflict
%        Oct. 17, 2008 - added VINFO support
%        Nov. 16, 2008 - history fix, doc update, code cleaning, rename
%                        from GENUMDESC to GETENUMDESC
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 16, 2008 at 05:15 GMT

% todo:

% require at least two inputs
if(nargin<2)
    error('seizmo:getenumdesc:notEnoughInputs',...
        'Not enough input arguments.')
end

% preallocate output
varnargin=length(varargin);
nvarargout=cell(1,varnargin);
varargout=nvarargout;
[varargout{:}]=deal(cell(numel(data),1));

% get header info
[nvarargout{:}]=getheader(data,varargin{:});

% loop over versions
[h,idx]=versioninfo(data);
for i=1:numel(h)
    % indexing of data with this header version
    ind=find(idx==i);
    
    % loop over fields
    for j=1:length(varargin)
        % check for cell output (char field)
        if(iscell(nvarargout{j}(ind)))
            error('seizmo:getenumdesc:badField',...
                'String fields are not supported!');
        end
        
        % compare
        bad=isnan(nvarargout{j}(ind));
        undef=nvarargout{j}(ind)==h(i).undef.ntype;
        unknown=~(bad | undef) & (nvarargout{j}(ind)<h(i).enum(1).minval...
                | nvarargout{j}(ind)>h(i).enum(1).maxval ...
                | fix(nvarargout{j}(ind))~=nvarargout{j}(ind));
        good=~(bad | undef | unknown);
        
        % assign enum descriptions
        if(any(bad))
            varargout{j}(ind(bad))={'Unknown Enum Field'};
        end
        if(any(undef))
            varargout{j}(ind(undef))={'Undefined Enum Field'};
        end
        if(any(unknown))
            varargout{j}(ind(unknown))={'Unknown Enum Value'};
        end
        if(any(good))
            varargout{j}(ind(good))=...
                h(i).enum(1).desc(nvarargout{j}(ind(good))+1);
        end
    end
end

end
