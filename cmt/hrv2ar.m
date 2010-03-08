function [varargout]=hrv2ar(varargin)
%HRV2AR    Convert moment tensor from Harvard from to Aki & Richards form
%
%    Usage:    momten=hrv2ar(momten)
%              [Mxx,Myy,Mzz,Mxy,Mxz,Myz]=hrv2ar(Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)
%
%    Description: MOMTEN=HRV2AR(MOMTEN) converts moment tensors stored in
%     MOMTEN from Harvard/USGS form (Up, South, East) to Aki & Richards
%     form (North, East, Up).  MOMTEN may be a 3x3xN array or a Nx6 array
%     where the N allows for multiple moment tensors.  The Nx6 array case
%     must have the layout [Mrr Mtt Mpp Mrt Mrp Mtp] which allows for
%     compact moment tensor storage compared to the 3x3xN case which
%     repeats Mrt, Mrp, Mtp elements in the upper triange.  The output
%     moment tensor will have the same size as that of the input.
%
%     [Mxx,Myy,Mzz,Mxy,Mxz,Myz]=HRV2AR(Mrr,Mtt,Mpp,Mrt,Mrp,Mtp) lets the
%     moment tensor elements be input and output separately.  There must be
%     exactly six inputs and they must be column vectors!
%
%    Notes:
%
%    Examples:
%     Converting to strike, dip, rake requires conversion of Harvard moment
%     tensors to Aki & Richard form:
%      sdr=mt2sdr(hrv2ar(mt));
%
%    See also: AR2HRV

%     Version History:
%        Mar.  8, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  8, 2010 at 13:50 GMT

% todo:

%AR  XX  YY  ZZ  XY  XZ  YZ
%     N   E   U  NE  NU  EU
%     2   3   1  -6   4  -5 (hrv2ar)
%
%HRV RR  TT  PP  RT  RP  TP
%     U   S   E  US  UE  SE
%     3   1   2   5  -6  -4 (ar2hrv)

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
    elseif(isequal(6,sz(2)))
        varargout{1}=varargin{1}(:,[2 3 1 6 4 5],:);
        varargout{1}(:,[4 6],:)=-varargout{1}(:,[4 6],:);
    else
        error('seizmo:hrv2ar:badInput',...
            'MOMTEN inproper size!');
    end
elseif(nargin==6)
    if(any(~cellfun('isreal',varargin) | cellfun('size',varargin,2)~=1))
        error('seizmo:hrv2ar:badInput',...
            'All inputs must be real-valued Nx1 vectors!');
    end
    [varargout{:}]=deal(varargin{[2 3 1 6 4 5]});
    varargout{4}=-varargout{4};
    varargout{6}=-varargout{6};
else
    error('seizmo:hrv2ar:badNumInput',...
        'Incorrect number of inputs!');
end

end
