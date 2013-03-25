function [report]=mt_check(varargin)
%MT_CHECK    Checks moment tensor data (array or struct)
%
%    Usage:    msg=mt_check(mt)
%              msg=mt_check(mrr,mtt,mpp,mrt,mrp,mtp)
%              msg=mt_check(mxx,myy,mzz,mxy,mxz,myz)
%
%    Description:
%     MSG=MT_CHECK(MT) checks that MT is a valid struct or array for moment
%     tensor functions.  The struct must be of the form output by FINDCMT
%     or FINDCMTS or a numeric array as:
%        - 1x6  [Mrr Mtt Mpp Mrt Mrp Mtp]
%        - 3x3  [Mrr Mrt Mrp; Mrt Mtt Mtp; Mrp Mtp Mpp]
%     assuming the components are Harvard style.  The numeric arrays must
%     be real-valued and may be Nx6 or 3x3xN to allow for multiple moment
%     tensor processing.
%
%     MSG=MT_CHECK(MRR,MTT,MPP,MRT,MRP,MTP)
%     MSG=MT_CHECK(MXX,MYY,MZZ,MXY,MXZ,MYZ) allow the individual components
%     to be checked.  This is generally not a useful moment tensor format.
%
%    Notes:
%
%    Examples:
%     % Check that FINDCMTS output is okay:
%     error(mt_check(findcmts))
%
%    See also: MT_CHANGE, FINDCMT, FINDCMTS

%     Version History:
%        Mar. 24, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2013 at 13:30 GMT

% todo:

% check nargin
report=[];
error(nargchk(1,6,nargin));

% act by number of inputs
if(nargin==1) % scalar struct or nx6/3x3xn array
    mt=varargin{1};
    VALID_CMT_FIELDS={'scalarmoment' 'exponent' 'year' 'month' 'day' ...
        'hour' 'minute' 'seconds' 'centroidtime' 'centroidlat' 'name' ...
        'centroidlon' 'centroiddep' 'latitude' 'longitude' 'depth' ...
        'mrr' 'mtt' 'mpp' 'mrt' 'mrp' 'mtp'};
    if(isstruct(mt) && isscalar(mt))
        if(~all(ismember(VALID_CMT_FIELDS,fieldnames(mt))))
            report.identifier='seizmo:mt_check:badInput';
            report.message=...
                ['MT struct must have the following fields:\n' ...
                sprintf('%s ',VALID_CMT_FIELDS{:})];
            return;
        end
    elseif(isnumeric(mt) && isreal(mt))
        sz=size(mt);
        if(~(numel(sz)==2 && sz(2)==6) ...
                && ~(any(numel(sz)==[2 3]) && all(sz(1:2)==3)))
            report.identifier='seizmo:mt_check:badInput';
            report.message='MT array must be Nx6 or 3x3xN!';
            return;
        end
    else
        report.identifier='seizmo:mt_check:badInput';
        report.message=['MT must be a scalar struct or ' ...
            'a Nx6 / 3x3xN real array!'];
        return;
    end
elseif(nargin==6) % components
    [mrr,mtt,mpp,mrt,mrp,mtp]=deal(varargin{:});
    if(~isnumeric(mrr) || ~isreal(mrr) || ~isnumeric(mtt) ...
            || ~isreal(mtt) || ~isnumeric(mpp) || ~isreal(mpp) ...
            || ~isnumeric(mrt) || ~isreal(mrt) || ~isnumeric(mrp) ...
            || ~isreal(mrp) || ~isnumeric(mtp) || ~isreal(mtp))
        report.identifier='seizmo:mt_check:badInput';
        report.message=['Moment tensor components must be numeric ' ...
            'and real-valued!'];
        return;
    elseif(~isequalsizeorscalar(varargin{:}))
        report.identifier='seizmo:mt_check:badInput';
        report.message=['Moment tensor components must be scalar or ' ...
            'equally sized!'];
        return;
    end
else
    report.identifier='seizmo:mt_check:badInput';
    report.message=['Moment tensor(s) must be specified ' ...
        'with 1 or 6 inputs!'];
    return;
end

end
