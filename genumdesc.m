function [varargout]=genumdesc(data,varargin)
%GENUMDESC    Get enum string description from SAClab header enum field
%
%    Description: GENUMDESC(DATA,FIELD) returns a cellstring array filled 
%     with the description strings associated with the enum field FIELD
%     stored in the SAClab structure DATA.  This is more useful/readible
%     than the magic number (enumerator integer) returned with GH.
%
%     GENUMDESC(DATA,FIELD1,...,FIELDN) returns a cellstring array for each
%     field supplied.
%
%    Notes:
%     - Numeric header fields not defined to be enum fields can be used as 
%       if they were enum fields with GENUMDESC.  This gives the user the 
%       ability to have more enum fields if needed.  Character fields are
%       NOT able to be treated as enum fields.
%     - Nonexistent header fields will return 'Unknown Enum Field'.
%     - Enum values that aren't defined in the enumerator list are tagged
%       as 'Unknown Enum Value' unless they are the defined UNDEFINED
%       value in which case they are tagged as 'Undefined Enum Field'.
%
%    Usage: cellstr=genumdesc(data,'field')
%           [cellstr1,...,cellstrN]=genumdesc(data,'field1',...,'fieldN')
%
%    Examples:
%     To check if all records are timeseries data:
%      if(all(strcmp(genumdesc(data,'iftype'),'Time Series File')))
%          disp('timeseries dataset')
%      end
%
%     Interpret resp0 as an enum field:
%      my_enum_desc=genumdesc(data,'resp0')   
%
%    See also: gh, glgc, genum

%     Version History:
%        ????????????? - Initial Version
%        June 13, 2008 - Documentation update, compat fixes, better
%                        support for undefined and unknown values/fields
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 13, 2008 at 20:15 GMT

% require at least two inputs
if(nargin<2)
    error('MATLAB:nargchk:notEnoughInputs',...
        'Not enough input arguments.')
end

% preallocate output
varnargin=length(varargin);
nvarargout=cell(1,varnargin);
varargout=nvarargout;
[varargout{:}]=deal(cell(length(data),1));

% get header info
[nvarargout{:}]=gh(data,varargin{:});

% loop over versions
v=[data.version];
for i=unique(v)
    % grab header setup
    h=seisdef(i);
    
    % indexing of data with this header version
    ind=find(v==i);
    
    % loop over fields
    for j=1:length(varargin)
        % unknown enum field if NaN
        nan=isnan(nvarargout{j}(ind));
        if(any(nan))
            varargout{j}(ind(nan))={'Unknown Enum Field'};
        end
        
        % undefined enum value
        undef=nvarargout{j}(ind)==h.undef.ntype;
        if(any(undef))
            varargout{j}(ind(undef))={'Undefined Enum Field'};
        end
        
        % unknown enum value if outside enumerator list
        unknown=(~nan & ~undef & (nvarargout{j}(ind)<h.enum(1).minval ...
                | nvarargout{j}(ind)>h.enum(1).maxval ...
                | round(nvarargout{j}(ind))~=nvarargout{j}(ind)));
        if(any(unknown))
            varargout{j}(ind(unknown))={'Unknown Enum Value'};
        end
        
        % good enum values
        good=~nan & ~undef & ~unknown;
        if(any(good))
            varargout{j}(ind(good))=...
                h.enum(1).desc(nvarargout{j}(ind(good))+1);
        end
    end
end

end
