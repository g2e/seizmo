function [varargout]=getenumid(data,varargin)
%GETENUMID    Get enum string id from SEIZMO data header enum field
%
%    Usage: cellstr=getenumid(data,'field')
%           [cellstr1,...,cellstrN]=getenumid(data,'field1',...,'fieldN')
%
%    Description: GETENUMID(DATA,FIELD) returns a cellstring array filled
%     with the enum id strings associated with the enum field FIELD stored
%     in the SEIZMO structure DATA.  This is more useful/readible than the 
%     magic number (enumerator integer) returned with GETHEADER.
%
%     GETENUMID(DATA,FIELD1,...,FIELDN) returns a cellstring array for each
%     field supplied.
%
%    Notes:
%     - Numeric header fields not defined to be enum fields can be used as 
%       if they were enum fields with GETENUMID.  This gives the user the 
%       ability to have more enum fields if needed.  Character fields are
%       NOT able to be treated as enum fields.
%     - Nonexistent header fields will return 'Unknown Enum Field'.
%     - Enum values that aren't defined in the enumerator list are tagged
%       as 'Unknown Enum Value' unless they are the defined UNDEFINED
%       value in which case they are tagged as 'Undefined Enum Field'. This
%       avoids conflict with enum fields called 'unknown' and 'undefined'.
%
%    Examples:
%     To check if all records are timeseries data:
%      if(all(strcmp(getenumid(data,'iftype'),'itime')))
%          disp('timeseries dataset')
%      end
%
%     Interpret resp0 as an enum field:
%      my_enum_id=getenumid(data,'resp0')   
%
%    See also: getheader, getlgc, getenumdesc

%     Version History:
%        Feb. 23, 2008 - initial version
%        Feb. 28, 2008 - minor code clean
%        Feb. 29, 2008 - handle undefined
%        Mar.  4, 2008 - minor doc update
%        June 13, 2008 - doc update, compat fixes, handle values out of
%                        range or non-whole, avoid name conflict
%        Oct. 17, 2008 - added VINFO support
%        Nov. 16, 2008 - history fix, doc update, code cleaning, rename
%                        from GENUM to GETENUMID
%        Apr. 23, 2009 - move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:45 GMT

% require at least two inputs
if(nargin<2)
    error('seizmo:getenumid:notEnoughInputs',...
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
            error('seizmo:getenumid:badField',...
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
                h(i).enum(1).id(nvarargout{j}(ind(good))+1);
        end
    end
end

end
