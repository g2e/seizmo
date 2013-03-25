function [varargout]=ar2hrv(varargin)
%AR2HRV    Convert moment tensor from Aki & Richards form to Harvard form
%
%    Usage:    mt=ar2hrv(mt)
%              mt=ar2hrv(Mxx,Myy,Mzz,Mxy,Mxz,Myz)
%              [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp]=ar2hrv(...)
%
%    Description:
%     MT=AR2HRV(MT) converts moment tensors stored in MT from Aki &
%     Richards form (North, East, Down) to Harvard/USGS form (Up, South,
%     East).  MT may be a 3x3xN array arranged as:
%        [Mxx Mxy Mxz
%         Mxy Myy Myz
%         Mxz Myz Mzz]
%     or a Nx6 array as:
%        [Mxx Myy Mzz Mxy Mxz Myz]
%     where N allows for multiple moment tensors down the corresponding
%     dimension.  The output moment tensor will have the same size as that
%     of the input.
%
%     MT=AR2HRV(Mxx,Myy,Mzz,Mxy,Mxz,Myz) allows inputing the components
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
%        Mar. 25, 2013 - update for mt_change/mt_check
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 13:50 GMT

% todo:

%HRV RR  TT  PP  RT  RP  TP
%    UU  SS  EE  US  UE  SE
%     3   1   2   5  -6  -4 (ar2hrv)
%AR  XX  YY  ZZ  XY  XZ  YZ
%    NN  EE  DD  NE  ND  ED
%     2   3   1  -6   4  -5 (hrv2ar)

% nargin check
error(nargchk(1,6,nargin));

% check moment tensor and force to vector form
error(mt_check(varargin{:}));
[mt,from]=mt_change('v',varargin{:});
if(from=='s')
    error('seizmo:ar2hrv:badInput',...
        'AR2HRV does not allow scalar struct MT input!');
end

% ar2hrv
mt=mt(:,[3 1 2 5 6 4]);
mt(:,[5 6])=-mt(:,[5 6]);

% output/separate
if(nargout>1)
    varargout=cell(1,6);
    [varargout{:}]=mt_change('c',mt);
elseif(from=='g')
    varargout={mt_change('g',mt)};
else
    varargout={mt};
end

end
