function [varargout]=genumdesc(data,varargin)
%GENUMDESC    Get enum string description from SAClab enum header field
%
%    Description: Returns cellstrings containing enum description strings 
%     corresponding to each field/value.
%
%    Usage: [enumcellstr1,enumcellstr2,...]=genum(data,'field1','field2',...)
%
%    Examples:
%     To check if all records are timeseries:
%      if(all(strcmp(genumdesc(data,'iftype'),'Time Series File')))
%          disp('timeseries dataset')
%      end
%
%    See also: gh, glgc, genum

% do nothing on no input
if(nargin<2); return; end

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
    h=sachi(i);
    
    % indexing of data with this header version
    ind=find(v==i);
    
    % loop over fields
    for j=1:length(varargin)
        varargout{j}(ind)=h.enum(1).desc(nvarargout{j}(ind)+1);
    end
end

end

