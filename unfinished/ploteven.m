function [varargout]=ploteven(data,varargin)
%PLOTEVEN    Plots SEIZMO records evenly spaced in the vertical direction
%
%    Usage:
%
%    Description:
%
%    Notes:
%
%    Examples:
%
%    See also:

%     Version History:
%        Mar. 18, 2008 - complete rewrite, renamed from PLOT0 to PLOTEVEN
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 18, 2010 at 22:25 GMT

% todo:

% check data structure
[h,idx]=versioninfo(data,'dep');

% parse options
opt=plotparameters(false,varargin{:});

% verbosity
verbose=seizmoverbose;

% number of records
nrecs=numel(data);

% select/open figure & axis
if(~isempty(opt.SPECIAL.HANDLE) && ishandle(opt.SPECIAL.HANDLE))
    if(opt.SPECIAL.HANDLE==fix(opt.SPECIAL.HANDLE))
        % select figure
        varargout{1}=figure(opt.SPECIAL.HANDLE);
    else
        % select axis
        varargout{1}=get(opt.SPECIAL.HANDLE,'parent');
        subplot(opt.SPECIAL.HANDLE);
    end
else
    % open new figure
    varargout{1}=figure;
end

% style the plot
set(gcf,opt.FIGURE{:});
set(gca,opt.AXES{:});

% record coloring
try
    % try using as a colormap function
    cmap=str2func(opt.SPECIAL.COLORMAP);
    colors=cmap(nrecs);
catch
    % otherwise it is a color array/string (each row gives a color)
    % - repeat vertically until it exceeds nrecs
    colors=repmat(opt.SPECIAL.COLORMAP,...
        ceil(nrecs/size(opt.SPECIAL.COLORMAP,1)),1);
end

% get header info
even=~strcmpi(getlgc(data,'leven'),'false');
iftype=getenumid(data,'iftype');
[t,kt,o,ko,a,ka,f,kf,b,e,npts,delta,depmin,depmax]=getheader(data,...
    't','kt','o','ko','a','ka','f','kf','b','e','npts','delta',...
    'depmin','depmax');

% extract the desired component for spectral records


end
