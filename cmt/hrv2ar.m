function [varargout]=hrv2ar(varargin)
%HRV2AR    Convert moment tensor from Harvard from to Aki & Richards form
%
%    Usage:    momten=hrv2ar(momten)
%              momten=hrv2ar(Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)
%              [Mxx,Myy,Mzz,Mxy,Mxz,Myz]=hrv2ar(...)
%
%    Description:
%     MOMTEN=HRV2AR(MOMTEN) converts moment tensors stored in MOMTEN from
%     Harvard/USGS form (Up, South, East) to Aki & Richards form (North,
%     East, Down).  MOMTEN may be a 3x3xN array arranged as:
%        [Mrr Mrt Mrp
%         Mrt Mtt Mtp
%         Mrp Mtp Mpp]
%     or a Nx6 array as:
%        [Mrr Mtt Mpp Mrt Mrp Mtp]
%     where N allows for multiple moment tensors down the corresponding
%     dimension.  The output moment tensor will have the same size as that
%     of the input.
%
%     MOMTEN=HRV2AR(Mrr,Mtt,Mpp,Mrt,Mrp,Mtp) allows inputing the components
%     of the moment tensor separately.  Output is an Nx6 array.  There must
%     be exactly six inputs and they must be column vectors or scalars!
%
%     [Mxx,Myy,Mzz,Mxy,Mxz,Myz]=HRV2AR(...) outputs the moment tensor
%     components separately as 6 Nx1 column vectors.
%
%    Notes:
%
%    Examples:
%     % Note the scalar moment and moment magnitude does not change:
%     [momentmag(mt) momentmag(hrv2ar(mt))]
%     [scalarmoment(mt) scalarmoment(hrv2ar(mt))]
%
%    See also: AR2HRV

%     Version History:
%        Mar.  8, 2010 - initial version
%        June  1, 2011 - doc update, improved usage
%        Mar. 13, 2013 - major doc fixes (A&R is NED not NEU!)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2013 at 13:50 GMT

% todo:

%HRV RR  TT  PP  RT  RP  TP
%    UU  SS  EE  US  UE  SE
%     3   1   2   5  -6  -4 (ar2hrv)
%AR  XX  YY  ZZ  XY  XZ  YZ
%    NN  EE  DD  NE  ND  ED
%     2   3   1  -6   4  -5 (hrv2ar)

if(nargin==1)
    if(~isreal(varargin{1}))
        error('seizmo:hrv2ar:badInput',...
            'MOMTEN must be real-valued!');
    end
    sz=size(varargin{1});
    if(isequal([3 3],sz(1:2)))
        varargout{1}=varargin{1}([2 3 1],[2 3 1],:);
        varargout{1}([1 3],:,:)=-varargout{1}([1 3],:,:);
        varargout{1}(:,[1 3],:)=-varargout{1}(:,[1 3],:);
        gflag=true;
    elseif(isequal(6,sz(2)))
        varargout{1}=varargin{1}(:,[2 3 1 6 4 5]);
        varargout{1}(:,[4 6])=-varargout{1}(:,[4 6]);
        gflag=false;
    else
        error('seizmo:hrv2ar:badInput',...
            'MOMTEN inproper size!');
    end
elseif(nargin==6)
    if(any(~cellfun('isreal',varargin) | cellfun('size',varargin,2)~=1))
        error('seizmo:hrv2ar:badInput',...
            'All inputs must be real-valued Nx1 column vectors!');
    end
    % expand scalars
    n=cellfun('prodofsize',varargin);
    sz=size(varargin{find(n==max(n),1,'first')});
    for i=1:6
        if(n(i)==1)
            varargin{i}=varargin{i}(ones(sz));
        else
            if(~isequal(size(varargin{i}),sz))
                error('seizmo:hrv2ar:badInput',...
                    'Non-scalar inputs must be equal sized!');
            end
        end
    end
    varargout{1}=[varargin{[2 3 1 6 4 5]}];
    varargout{1}(:,[4 6])=-varargout{1}(:,[4 6]);
    gflag=false;
else
    error('seizmo:hrv2ar:badNumInput',...
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
