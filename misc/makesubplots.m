function [ax]=makesubplots(r,c,p,varargin)
%MAKESUBPLOTS    Makes subplots in current figure returning axes handles
%
%    Usage:    ax=makesubplots(r,c)
%              ax=makesubplots(r,c,p)
%              ax=makesubplots(r,c,p,...)
%
%    Description:
%     AX=MAKESUBPLOTS(R,C) initializes all subplots in a RxC array of axes
%     in the current figure.  R specifies the number of rows and C gives
%     the number of columns.  So AX will be a CxR matrix of axes handles
%     (note this is actually transposed from the way the plots are arranged
%     -- so AX(2) corresponds to subplot(R,C,2) etc).  Transpose AX to get
%     axes handles in the visual ordering.  If no figure exists, one is
%     initialized.
%
%     AX=MAKESUBPLOTS(R,C,P) initializes only the subplots at positions P.
%     Note that this does not support subplots spanning multiple rows
%     and/or columns.  AX is a column vector in this case.
%
%     AX=MAKESUBPLOTS(R,C,P,...) allows passing additional options to
%     SUBPLOT.
%
%    Notes:
%
%    Examples:
%     % make a figure with 4 2x2 groups of subplots and add super
%     % labeling and super colorbars to each group
%     fh=figure;
%     set(fh,'position',get(fh,'position').*[1 1 1.5 1.5]);
%     ax=makesubplots(5,5,submat(lind(5),1:2,[1 2 4 5]),'parent',fh);
%     ax=mat2cell(reshape(ax,4,4),[2 2],[2 2]);
%     for i=1:4
%         supertitle(ax{i},['super title ' num2str(i)]);
%         superxlabel(ax{i},['super xlabel ' num2str(i)]);
%         superylabel(ax{i},['super ylabel ' num2str(i)]);
%         supercolorbar(ax{i},'location','eastoutside');
%     end
%
%    See also: SUBPLOT, FIGURE, AXEXPAND, AXMOVE, AXSTRETCH, NOLABELS,
%              NOTICKS, NOTITLES, NOCOLORBARS, SUPERTITLE, SUPERXLABEL,
%              SUPERYLABEL, SUPERCOLORBAR

%     Version History:
%        Aug.  4, 2010 - initial version
%        Aug.  8, 2010 - bugfix: support p as matrix, logical support
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 12:25 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));

% check rows/cols
if(~isreal(r) || ~isscalar(r) || r~=fix(r) || r<1)
    error('seizmo:makesubplots:badInput',...
        'R must be a positve integer!');
elseif(~isreal(c) || ~isscalar(c) || c~=fix(c) || c<1)
    error('seizmo:makesubplots:badInput',...
        'C must be a positve integer!');
end

% default p
axarray=false;
if(nargin<3 || isempty(p)); p=1:r*c; axarray=true; end

% check p
if(islogical(p))
    if(~isequal(size(p),[r c]))
        error('seizmo:makesubplots:badInput',...
            'P must be an array of indices (logical or linear)!');
    end
    p=find(p');
elseif(~isreal(p) || any(p(:)~=fix(p(:))) || any(p(:)<1 | p(:)>r*c))
    error('seizmo:makesubplots:badInput',...
        'P must be an array of indices (logical or linear)!');
end

% setup output
np=numel(p);
if(axarray)
    ax=nan(c,r);
else
    ax=nan(np,1);
end

% loop over p
for i=1:np
    ax(i)=subplot(r,c,p(i),varargin{:});
end

end
