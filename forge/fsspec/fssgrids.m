function [baz,smag,es,ns]=fssgrids(s)
%FSSGRIDS    Returns the backazimuth and slowness grids for fss spectra
%
%    Usage:    [baz,smag,es,ns]=fssgrids(s)
%
%    Description:
%     [BAZ,SMAG,ES,NS]=FSSGRIDS(S) returns the backazimuth and slowness
%     grids for the frequency-slowness spectra in fss struct S.  S must
%     correspond to the output from function FSS.  BAZ, SMAG, ES, & NS are
%     all equal sized 2D arrays equal in size to the first 2 dimensions of
%     the spectra in S.
%
%    Notes:
%     - If S has multiple elements the outputs are cell arrays of
%       equal size to S where each cell gives the corresponding
%       grid.
%
%    Examples:
%     % Apply a gaussian filter based on slowness magnitude:
%     [baz,smag]=fssgrids(s);
%     x=gaussiantf(smag,30,5);
%     s.spectra=s.spectra.*repmat(x,[1 1 size(s.spectra,3)]);
%
%    See also: FSSIDX, FSSDBINFO, FSSSUB, FSSAVG,
%              FSS, FSSXC, FSSHORZ, FSSHORZXC

%     Version History:
%        Sep. 12, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 12, 2012 at 14:25 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check fk struct
error(chkfss(s));

% loop over every spectra
sz=size(s);
[baz,smag,es,ns]=deal(cell(sz));
for i=1:numel(s)
    % get grids of backazimuth and slowness
    if(s(i).polar)
        nx=numel(s(i).x);
        ny=numel(s(i).y);
        baz{i}=s(i).x(ones(ny,1),:);
        smag{i}=s(i).y(:,ones(nx,1));
        [es{i},ns{i}]=slowbaz2kxy(smag{i},baz{i});
    else % cartesian
        nx=numel(s(i).x);
        ny=numel(s(i).y);
        es{i}=s(i).x(ones(ny,1),:);
        ns{i}=s(i).y(:,ones(nx,1));
        [smag{i},baz{i}]=kxy2slowbaz(es{i},ns{i});
    end
end

% uncell scalar spectra
if(isscalar(s))
    baz=baz{1};
    smag=smag{1};
    es=es{1};
    ns=ns{1};
end

end
