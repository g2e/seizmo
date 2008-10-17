function [varargout]=genum(data,varargin)
%GENUM    Get enum string id from SAClab data header enum field
%
%    Description: GENUM(DATA,FIELD) returns a cellstring array filled with
%     the enum id strings associated with the enum field FIELD stored in 
%     the SAClab structure DATA.  This is more useful/readible than the 
%     magic number (enumerator integer) returned with GH.
%
%     GENUM(DATA,FIELD1,...,FIELDN) returns a cellstring array for each
%     field supplied.
%
%    Notes:
%     - Numeric header fields not defined to be enum fields can be used as 
%       if they were enum fields with GENUM.  This gives the user the 
%       ability to have more enum fields if needed.  Character fields are
%       NOT able to be treated as enum fields.
%     - Nonexistent header fields will return 'Unknown Enum Field'.
%     - Enum values that aren't defined in the enumerator list are tagged
%       as 'Unknown Enum Value' unless they are the defined UNDEFINED
%       value in which case they are tagged as 'Undefined Enum Field'.
%
%    Usage: cellstr=genum(data,'field')
%           [cellstr1,...,cellstrN]=genum(data,'field1',...,'fieldN')
%
%    Examples:
%     To check if all records are timeseries data:
%      if(all(strcmp(genum(data,'iftype'),'itime')))
%          disp('timeseries dataset')
%      end
%
%     Interpret resp0 as an enum field:
%      my_enum_id=genum(data,'resp0')   
%
%    See also: gh, glgc, genumdesc

%     Version History:
%        ????????????? - Initial Version
%        June 13, 2008 - Documentation update, compat fixes, better
%                        support for undefined and unknown values/fields
%        Oct. 17, 2008 - added VINFO support
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 17, 2008 at 02:45 GMT

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
[h,idx]=vinfo(data);
for i=1:numel(h)
    % indexing of data with this header version
    ind=find(idx==i);
    
    % loop over fields
    for j=1:length(varargin)
        % unknown enum field if NaN
        nan=isnan(nvarargout{j}(ind));
        if(any(nan))
            varargout{j}(ind(nan))={'Unknown Enum Field'};
        end
        
        % undefined enum value
        undef=nvarargout{j}(ind)==h(i).undef.ntype;
        if(any(undef))
            varargout{j}(ind(undef))={'Undefined Enum Field'};
        end
        
        % unknown enum value if outside enumerator list
        unknown=(~nan & ~undef & (nvarargout{j}(ind)<h(i).enum(1).minval...
                | nvarargout{j}(ind)>h(i).enum(1).maxval ...
                | round(nvarargout{j}(ind))~=nvarargout{j}(ind)));
        if(any(unknown))
            varargout{j}(ind(unknown))={'Unknown Enum Value'};
        end
        
        % good enum values
        good=~nan & ~undef & ~unknown;
        if(any(good))
            varargout{j}(ind(good))=...
                h(i).enum(1).id(nvarargout{j}(ind(good))+1);
        end
    end
end

end
