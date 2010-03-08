function [varargout]=ar2hrv(varargin)
%AR2HRV    Convert moment tensor from Aki & Richards form to Harvard form
%
%    Usage:    momten=ar2hrv(momten)
%              [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp]=ar2hrv(Mxx,Myy,Mzz,Mxy,Mxz,Myz)
%
%    Description: MOMTEN=AR2HRV(MOMTEN) converts moment tensors stored in
%     MOMTEN from Aki & Richards form (North, East, Up) to Harvard/USGS
%     form (Up, South, East).  MOMTEN may be a 3x3xN array or a Nx6 array
%     where the N allows for multiple moment tensors.  The Nx6 array case
%     must have the layout [Mxx Myy Mzz Mxy Mxz Myz] which allows for
%     compact moment tensor storage compared to the 3x3xN case which
%     repeats Mxy, Mxz, Myz elements in the upper triange.  The output
%     moment tensor will have the same size as that of the input.
%     
%     [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp]=AR2HRV(Mxx,Myy,Mzz,Mxy,Mxz,Myz) lets the
%     moment tensor elements be input and output separately.  There must be
%     exactly six inputs and they must be column vectors!
%
%    Notes:
%
%    Examples:
%     Plotting the elementary moment tensors properly requires conversion
%     to Harvard form (if you are using the BB function):
%      mt=elementary_mt(1:6);
%      figure; bb(ar2hrv(mt),1:6,ones(6,1),0.5*ones(6,1),0,'b');
%      axis equal off
%
%    See also: HRV2AR

%     Version History:
%        Mar.  8, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  8, 2010 at 13:50 GMT

% todo:

%AR  XX  YY  ZZ  XY  XZ  YZ
%     N   E   U  NE  NU  EU
%     2   3   1   6  -4  -5 (hrv2ar)
%
%HRV RR  TT  PP  RT  RP  TP
%     U   S   E  US  UE  SE
%     3   1   2   5  -6  -4 (ar2hrv)

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
    elseif(isequal(6,sz(2)))
        varargout{1}=varargin{1}(:,[3 1 2 5 6 4],:);
        varargout{1}(:,[5 6],:)=-varargout{1}(:,[5 6],:);
    else
        error('seizmo:ar2hrv:badInput',...
            'MOMTEN inproper size!');
    end
elseif(nargin==6)
    if(any(~cellfun('isreal',varargin) | cellfun('size',varargin,2)~=1))
        error('seizmo:ar2hrv:badInput',...
            'All inputs must be real-valued Nx1 vectors!');
    end
    [varargout{:}]=deal(varargin{[3 1 2 5 6 4]});
    varargout{5}=-varargout{5};
    varargout{6}=-varargout{6};
else
    error('seizmo:ar2hrv:badNumInput',...
        'Incorrect number of inputs!');
end

end
