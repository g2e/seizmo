function [ax]=makesubplots(r,c,p,varargin)
%MAKESUBPLOTS    Makes subplots in current figure returning axes handles
%
%    Usage:    ax=makesubplots(r,c)
%              ax=makesubplots(r,c,p)
%              ax=makesubplots(r,c,p,...)
%
%    Description: AX=MAKESUBPLOTS(R,C) initializes all subplots in a RxC
%     array of axes in the current figure.  R specifies the number of rows
%     and C gives the number of columns.  So AX will be a CxR matrix of
%     axes handles (note this is actually transposed from the way the plots
%     are arranged -- so AX(2) corresponds to subplot(R,C,2) etc).
%     Transpose AX to get axes handles in the visual ordering.  If no
%     figure exists, one is initialized.
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
%     % make a figure with 4x3 arrangement of subplots, then expand them by
%     % 15% and drop labels on any axes not at the figure edge
%     figure;
%     ax=makesubplots(5,3,1:12);
%     ax=reshape(ax,3,4);
%     tax=ax';
%     axexpand(ax,15);
%     nolabels(tax(5:12),'y');
%     nolabels(ax(1:9),'x');
%     noticks(tax(5:12),'y');
%     noticks(ax(1:9),'x');
%     th=supertitle(ax,'This is a sooooooooooooooooooooooooper title!');
%     ax0=get(th,'parent'); % make title,colorbar,ylabel share same axis
%     superxlabel(ax0,'This is a sooooooooooooooooooooooooper xlabel!');
%     superylabel(ax0,'This is a sooooooooooooooooooooooooper ylabel!');
%     cb=supercolorbar(ax,'location','south');
%     cpos=get(cb,'position');
%     set(cb,'position',[cpos(1) cpos(2)-.15 cpos(3) cpos(4)/2]);
%     set(cb,'xaxislocation','bottom');
%
%    See also: SUBPLOT, FIGURE, AXEXPAND, AXMOVE, NOLABELS, NOTICKS,
%              SUPERTITLE, SUPERXLABEL, SUPERYLABEL, SUPERCOLORBAR,
%              COMPACTAXES

%     Version History:
%        Aug.  4, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  4, 2010 at 12:25 GMT

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
if(~isreal(p) || any(p~=fix(p)) || any(p<1 | p>r*c))
    error('seizmo:makesubplots:badInput',...
        'P must be an array of positive integers!');
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
