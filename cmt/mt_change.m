function [varargout]=mt_change(to,varargin)
%MT_CHANGE    Alter moment tensor format (array or struct
%
%    Usage:    [mt,from]=mt_change(to,mt)
%              [mt,from]=mt_change(to,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)
%              [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp,from]=mt_change('c',mt)
%
%    Description:
%     [MT,FROM]=MT_CHANGE(TO,MT) allows converting from a scalar struct as
%     output by FINDCMT/FINDCMTS into a 3x3xN or Nx6 array or between the
%     3x3xN and Nx6 array formats.  Converting to a scalar struct is not
%     allowed as there is far more information in that format than in the
%     others.  TO must be either 'g' (returns 3x3xN array) or 'v' (returns
%     a Nx6 array).  FROM is either 's' (struct input), 'g' (3x3xN input),
%     or 'v' (Nx6 input).
%
%     [MT,FROM]=MT_CHANGE(TO,MRR,MTT,MPP,MRT,MRP,MTP) allows converting
%     individual moment tensor components into either a 3x3xN or Nx6 array.
%     FROM will be 'c' indicating that MT was created from individual
%     component inputs.
%
%     [MRR,MTT,MPP,MRT,MRP,MTP,FROM]=MT_CHANGE('C',MT) returns the moment
%     tensor components individually by specifying the first input as 'C'.
%
%    Notes:
%
%    Examples:
%     % Convert FINDCMTS output into a compact moment tensor format:
%     mt=mt_change('v',findcmts);
%
%    See also: MT_CHECK, FINDCMT, FINDCMTS

%     Version History:
%        Mar. 25, 2013 - initial version, allow "conversion" to same type
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 13:30 GMT

% todo:

% check nargin
error(nargchk(2,7,nargin));

% check output string
if(~isscalar(to) || ~ischar(to))
    error('seizmo:mt_change:badInput',...
        'TO must be a single character (''v'' ''g'' or ''c'')!');
elseif(~any(to=='cvg'))
    error('seizmo:mt_change:badInput',...
        'Unrecognized TO character (must be ''v'' ''g'' or ''c'')!');
end

% check input mt
error(mt_check(varargin{:}));

% act by number of inputs
if(nargin==2) % scalar struct or nx6/3x3xn array
    % what is the input type
    if(isstruct(varargin{1})) % struct
        switch to
            case 'v'
                % Nx6
                varargout{1}=[varargin{1}.mrr varargin{1}.mtt ...
                    varargin{1}.mpp varargin{1}.mrt varargin{1}.mrp ...
                    varargin{1}.mtp];
            case 'g'
                % 3x3xN
                varargout{1}=permute(cat(3,...
                    [varargin{1}.mrr varargin{1}.mrt varargin{1}.mrp],...
                    [varargin{1}.mrt varargin{1}.mtt varargin{1}.mtp],...
                    [varargin{1}.mrp varargin{1}.mtp varargin{1}.mpp]),...
                    [2 3 1]);
            case 'c'
                varargout=cell(1,6);
                [varargout{1:6}]=deal(varargin{1}.mrr,varargin{1}.mtt,...
                    varargin{1}.mpp,varargin{1}.mrt,varargin{1}.mrp,...
                    varargin{1}.mtp);
        end
        varargout{end+1}='s';
    elseif(size(varargin{1},2)==6) % vector
        switch to
            case 'v'
                varargout=varargin;
            case 'g'
                varargout{1}=reshape(...
                    varargin{1}(:,[1 4 5 4 2 6 5 6 3])',...
                    [3 3 size(varargin{1},1)]);
            case 'c'
                varargout=num2cell(varargin{1},1);
        end
        varargout{end+1}='v';
    else % grid
        switch to
            case 'v'
                varargin{1}=permute(varargin{1},[3 1 2]);
                varargout{1}=[varargin{1}(:,1) varargin{1}(:,5) ...
                    varargin{1}(:,9) varargin{1}(:,2) ...
                    varargin{1}(:,3) varargin{1}(:,6)];
            case 'g'
                varargout=varargin;
            case 'c'
                varargin{1}=permute(varargin{1},[3 1 2]);
                varargout={varargin{1}(:,1) varargin{1}(:,5) ...
                    varargin{1}(:,9) varargin{1}(:,2) ...
                    varargin{1}(:,3) varargin{1}(:,6)};
        end
        varargout{end+1}='g';
    end
elseif(nargin==7) % components
    % expand scalars & turn into column vectors
    [varargin{:}]=expandscalars(varargin{:});
    [varargin{:}]=mat2vec(varargin{:});
    
    % convert
    switch to
        case 'v'
            varargout{1}=cat(2,varargin{:});
        case 'g'
            varargin=cat(2,varargin{:});
            varargout{1}=reshape(varargin(:,[1 4 5 4 2 6 5 6 3])',...
                    [3 3 size(varargin,1)]);
        case 'c'
            varargout=varargin;
    end
    varargout{end+1}='c';
end

end
