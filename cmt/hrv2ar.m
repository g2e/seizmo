function [varargout]=hrv2ar(varargin)
%HRV2AR    Convert moment tensor from Harvard from to Aki & Richards form
%
%    Usage:    mt=hrv2ar(mt)
%              mt=hrv2ar(Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)
%              [Mxx,Myy,Mzz,Mxy,Mxz,Myz]=hrv2ar(...)
%
%    Description:
%     MT=HRV2AR(MT) converts moment tensors stored in MT from Harvard/USGS
%     form (Up, South, East) to Aki & Richards form (North, East, Down).
%     MT may be a 3x3xN array arranged as:
%        [Mrr Mrt Mrp
%         Mrt Mtt Mtp
%         Mrp Mtp Mpp]
%     or a Nx6 array as:
%        [Mrr Mtt Mpp Mrt Mrp Mtp]
%     where N allows for multiple moment tensors down the corresponding
%     dimension.  The output moment tensor will have the same size as that
%     of the input.
%
%     MT=HRV2AR(Mrr,Mtt,Mpp,Mrt,Mrp,Mtp) allows specifying the moment
%     tensor components separately.  Output is an Nx6 array.  There must
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
    error('seizmo:hrv2ar:badInput',...
        'HRV2AR does not allow scalar struct MT input!');
end

% hrv2ar
mt=mt(:,[2 3 1 6 4 5]);
mt(:,[4 6])=-mt(:,[4 6]);

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
