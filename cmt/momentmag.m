function [mw]=momentmag(mt)
%MOMENTMAG    Returns the moment magnitude for a moment tensor
%
%    Usage:    mw=momentmag(mt)
%
%    Description:
%     Mw=MOMENTMAG(MT) calculates the moment magnitude Mw for the moment
%     tensor(s) given in MT.  MT must be a GlobalCMT struct (formatted as
%     returned by FINDCMTS), a Nx6 array or a 3x3xN array.  Note that for
%     the numeric arrays the tensors should include the exponent.
%
%    Notes:
%
%    Examples:
%     % Calculate the moment magnitudes of all CMTs in the GlobalCMT
%     % project catalog and make a histogram:
%     hist(momentmag(findcmts),100);
%
%    See also: FINDCMTS, FINDCMT, SCALARMOMENT

%     Version History:
%        Mar. 11, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 11, 2011 at 23:55 GMT

% todo:

% check inputs
sz=size(mt);
if(isstruct(mt) && isscalar(mt)) % global cmt struct
    % use scalarmoment field as that is more accurate than using
    % the truncated moment tensor values which are always biased
    % to lower magnitudes
    mw=(2/3).*(log10(mt.scalarmoment.*10.^mt.exponent)-16.1);
    return;
elseif(numel(sz)==2 && sz(2)==6) % (Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)
    mo=scalarmoment(mt);
elseif(isequal(sz(1:2),[3 3])) % [Mrr Mrt Mrp; Mrt Mtt Mtp; Mrp Mtp Mpp]
    mo=scalarmoment(mt);
else
    error('seizmo:momentmag:badInput',...
        'MT is not in a valid format!');
end

% moment magnitude
mw=(2/3).*(log10(mo)-16.1);

end
