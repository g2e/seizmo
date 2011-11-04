function [X]=submat_eval(X,varargin)
%SUBMAT_EVAL    Returns a submatrix using eval
%
%    Usage:    Y=submat_eval(X,DIM1,LIST1,DIM2,LIST2,...)
%
%    Description:
%     Y=SUBMAT_EVAL(X,DIM,LIST) creates a matrix Y that is the matrix X
%     reduced along dimension DIM to the indices in LIST.  If DIM is a list
%     of dimensions, LIST is used to reduce each dimension.  DIM and LIST
%     may be strings that can be evaluated such as '1:end-1'.
%
%     Y=SUBMAT_EVAL(X,DIM1,LIST1,DIM2,LIST2,...) allows for access to
%     multiple dimensions independently.  Later inputs take indexing
%     preference if a dimension is indexed more than once.
%
%    Notes:
%
%    Examples:
%     % Return a matrix reduced on the last dimension to the last 3
%     % elements:
%     x=submat_eval(x,'end','end-2:end')
%
%     % Remove the last elements along all dimensions:
%     x=submat_eval(x,':','1:end-1')
%
%     % These are equivalent:
%     x=repmat(x,[1 2 ones(1,ndims(x)-2)])
%     x=submat_eval(x,2,'[1:end 1:end]')
%
%    See also: SUBMAT, EVAL, COLON OPERATOR (:), END, REPMAT

%     Version History:
%        Nov. 12, 2008 - initial version
%        Apr. 23, 2009 - move usage up
%        Sep.  8, 2009 - fixed error message
%        Sep. 21, 2009 - updated one-liner, removed unnecessary brackets
%        June 24, 2010 - error in example fixed
%        Aug.  9, 2010 - DIM now takes strings
%        Nov.  1, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  1, 2011 at 01:05 GMT

% todo:

% CHECK VARARGIN
if(~mod(nargin,2))
    error('misc:submat_eval:badNumArgs','Unpaired DIM,LIST!');
end

% FIX ':' DIMENSION
NsD=ndims(X); % non-singleton dimensions
varargin(2.*find(strcmp(':',varargin(1:2:end)))-1)={1:NsD};
for i=1:2:nargin-2
    if(ischar(varargin{i}))
        % replace the word 'end' with final non-singleton dimension
        varargin{i}=regexprep(varargin{i},'\<end\>',num2str(NsD));
        varargin{i}=eval(varargin{i});
    end
end

% DEFAULT TO ENTIRE MATRIX AND EXPAND TO MAX INPUT DIMENSION
list(1:max([ndims(X) [varargin{1:2:end}]]))={':'};

% REDUCTION/REPLICATION OF DIMENSIONS
for i=1:2:nargin-2
    if(ischar(varargin{i+1}))
        [list{varargin{i}}]=deal(varargin{i+1});
    else
        [list{varargin{i}}]=deal(mat2str(varargin{i+1}));
    end
end

% FORM STRING FOR EVAL
list(1)=strcat('X(',list(1));               % ADD FRONT
list(1:end-1)=strcat(list(1:end-1),',');    % ADD COMMAS
list(end)=strcat(list(end),')');            % ADD REAR
string=strcat(list{:});                     % COMBINE
clear list varargin i                       % CLEAN UP

% EVAL
X=eval(string);

end
