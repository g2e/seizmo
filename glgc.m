function [varargout]=glgc(data,varargin)
%GLGC    Get logical string from SAClab logical header field
%
%    Description: GLGC(DATA,FIELD) returns a cellstring array containing
%     'true' 'false' 'undefined' or 'unknown' corresponding to values of
%     the logical field FIELD in the SAClab structure DATA.  This provides
%     a consistent logic framework than the output of GH for working with 
%     multiple datafile versions.
%
%     GLGC(DATA,FIELD1,...,FIELDN) returns a cellstring array for each
%     field supplied.
%    
%    Notes:
%     - Numeric header fields not defined to be logical can be used with 
%       GLGC as if they were logical fields.  This gives the user the 
%       ability to have more logical fields if needed.  Character fields 
%       are NOT able to be treated as logical fields.
%     - Nonexistent header fields will return 'Unknown Logic Field'.
%     - Values that are not in the logical definition list are marked as
%       'unknown' unless they are the UNDEFINED value in which case they 
%       are tagged as 'undefined'.
%
%    Usage: [lgccellstr1,...,lgccellstrN]=glgc(data,'field1',...,'fieldN')
%
%    Examples:
%     To check if all records are evenly spaced:
%      if(all(strcmp(glgc(data,'leven'),'true'))) 
%          disp('evenly spaced data')
%      end
%
%     Treat field RESP0 as a logical field:
%      my_lgc=glgc(data,'resp0')
%
%    See also: gh, genum, genumdesc

%     Version History:
%        ????????????? - initial version
%        June 12, 2008 - doc update
%        June 13, 2008 - sorted out undefined and unknown values/fields
%        Sep. 28, 2008 - output cleanup
%        Oct. 17, 2008 - added VINFO support
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 17, 2008 at 02:45 GMT

% todo:
% - history fix
% - vinfo
% - doc update

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
        if(any(nvarargout{j}(ind)==h(i).true))
            varargout{j}(ind(nvarargout{j}(ind)==h(i).true))={'true'};
        end
        if(any(nvarargout{j}(ind)==h(i).false))
            varargout{j}(ind(nvarargout{j}(ind)==h(i).false))={'false'};
        end
        if(any(nvarargout{j}(ind)==h(i).undef.ntype))
            varargout{j}(ind(nvarargout{j}(ind)==h(i).undef.ntype))={'undefined'};
        end
        if(any(isnan(nvarargout{j}(ind))))
            varargout{j}(ind(isnan(nvarargout{j}(ind))))={'Unknown Logic Field'};
        end
        if(any(nvarargout{j}(ind)~=h(i).true ...
            & nvarargout{j}(ind)~=h(i).false ...
            & nvarargout{j}(ind)~=h(i).undef.ntype ...
            & ~isnan(nvarargout{j}(ind))))
            [varargout{j}{ind(nvarargout{j}(ind)~=h(i).true ...
                & nvarargout{j}(ind)~=h(i).false ...
                & nvarargout{j}(ind)~=h(i).undef.ntype ...
                & ~isnan(nvarargout{j}(ind)))}]=deal('unknown');
        end
    end
end

end
