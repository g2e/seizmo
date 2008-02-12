function [varargout]=glgc(data,varargin)
%GLGC    Get universal truth values of SAClab logic header field
%
%    Description: Returns 116 ('t') if logic field is true, 
%                         102 ('f') if false and 
%                         117 ('u') otherwise (undefined/unknown).
%
%    Usage: logic=glgc(data,'leven')
%
%    Examples:
%     to check if a record is evenly spaced
%     if(116==glgc(data(1),'leven')); disp('evenly spaced'); end
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: gh

% do nothing on no input
if(nargin<1); return; end

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head'))
    error('data structure does not have proper fields')
end

% number of files
nrecs=length(data);

% headers setup
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% loop over fields
for i=1:length(varargin)
    % get header values
    varargout{i}=gh(data,varargin{i});
end

% loop over records
for i=1:nrecs
    % header logical index
    v=data(i).version==vers;
    
    % loop over fields
    for j=1:length(varargin)
        % check logic
        if(varargout{j}(i)==h(v).true)
            varargout{j}(i)=116; % 't'
        elseif(varargout{j}(i)==h(v).false)
            varargout{j}(i)=102; % 'f'
        else
            varargout{j}(i)=117; % 'u'
        end
    end
end

end