function [varargout]=ar2hrv(varargin)
%AR2HRV    Convert moment tensor from Aki & Richards form to Harvard form
%
%    Usage:    momten=ar2hrv(momten)
%              momten=ar2hrv(Mxx,Myy,Mzz,Mxy,Mxz,Myz)
%              [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp]=ar2hrv(...)
%
%    Description:
%     MOMTEN=AR2HRV(MOMTEN) converts moment tensors stored in MOMTEN from
%     Aki & Richards form (North, East, Down) to Harvard/USGS form (Up,
%     South, East).  MOMTEN may be a 3x3xN array arranged as:
%        [Mxx Mxy Mxz
%         Mxy Myy Myz
%         Mxz Myz Mzz]
%     or a Nx6 array as:
%        [Mxx Myy Mzz Mxy Mxz Myz]
%     where N allows for multiple moment tensors down the corresponding
%     dimension.  The output moment tensor will have the same size as that
%     of the input.
%
%     MOMTEN=AR2HRV(Mxx,Myy,Mzz,Mxy,Mxz,Myz) allows inputing the components
%     of the moment tensor separately.  Output is an Nx6 array.  There must
%     be exactly six inputs and they must be column vectors or scalars!
%     
%     [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp]=AR2HRV(...) outputs the moment tensor
%     components separately as 6 Nx1 column vectors.
%
%    Notes:
%
%    Examples:
%     % An explosion is expressed the same in both systems:
%     ar2hrv([1 1 1 0 0 0])
%
%    See also: HRV2AR

%     Version History:
%        Mar.  8, 2010 - initial version
%        June  1, 2011 - doc update, improved usage
%        Mar. 13, 2013 - major doc fixes (A&R is NED not NEU!)
%        Mar. 23, 2013 - changed example as the old one was wrong
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 23, 2013 at 13:50 GMT

% todo:

%HRV RR  TT  PP  RT  RP  TP
%    UU  SS  EE  US  UE  SE
%     3   1   2   5  -6  -4 (ar2hrv)
%AR  XX  YY  ZZ  XY  XZ  YZ
%    NN  EE  DD  NE  ND  ED
%     2   3   1  -6   4  -5 (hrv2ar)

if(nargin==1)
    if(~isreal(varargin{1}))
        error('seizmo:ar2hrv:badInput',...
            'MOMTEN must be real-valued!');
    end
    sz=size(varargin{1});
    if(isequal([3 3],sz(1:2)))
        varargout{1}=varargin{1}([3 1 2],[3 1 2],:);
        varargout{1}(3,:,:)=-varargout{1}(3,:,:);
        varargout{1}(:,3,:)=-varargout{1}(:,3,:);
        gflag=true;
    elseif(isequal(6,sz(2)))
        varargout{1}=varargin{1}(:,[3 1 2 5 6 4],:);
        varargout{1}(:,[5 6],:)=-varargout{1}(:,[5 6],:);
        gflag=false;
    else
        error('seizmo:ar2hrv:badInput',...
            'MOMTEN inproper size!');
    end
elseif(nargin==6)
    if(any(~cellfun('isreal',varargin) | cellfun('size',varargin,2)~=1))
        error('seizmo:ar2hrv:badInput',...
            'All inputs must be real-valued Nx1 vectors!');
    end
    % expand scalars
    n=cellfun('prodofsize',varargin);
    sz=size(varargin{find(n==max(n),1,'first')});
    for i=1:6
        if(n(i)==1)
            varargin{i}=varargin{i}(ones(sz));
        else
            if(~isequal(size(varargin{i}),sz))
                error('seizmo:ar2hrv:badInput',...
                    'Non-scalar inputs must be equal sized!');
            end
        end
    end
    varargout{1}=[varargin{[3 1 2 5 6 4]}];
    varargout{1}(:,[5 6])=-varargout{1}(:,[5 6]);
    gflag=false;
else
    error('seizmo:ar2hrv:badNumInput',...
        'Incorrect number of inputs!');
end

% separate if desired
if(nargout>1)
    if(gflag) % matrix
        [varargout{1:6}]=mt_g2c(varargout{1});
    else % vector
        [varargout{1:6}]=mt_v2c(varargout{1});
    end
end

end
