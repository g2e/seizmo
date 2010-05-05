function [c]=invertcolor(c,flag)
%INVERTCOLOR    Inverts colors given as rgb triplet or as short/long names
%
%    Usage:    ci=invertcolor(c)
%              ci=invertcolor(c,flag)
%
%    Description: CI=INVERTCOLOR(C) returns the "opposite" rgb color values
%     CI to those in C.  This is essentially CI=1-C but with special
%     handling of short/long color names (where it returns the opposing
%     short/long color name).  Note that "gray" [0.5 0.5 0.5] here is the
%     opposite of itself and so will not provide any contrast (see next
%     usage form for high contrast colors).
%
%     CI=INVERTCOLOR(C,FLAG) toggles high contrast mode.  FLAG set to TRUE
%     will find the color CI with the highest contrast to C using the
%     algorithm CI=1-ROUND(C).  This provides high contrast for grays
%     as all colors returned in high contrast mode are fully (de)saturated.
%     The default for FLAG is FALSE (return rgb color inverse).
%
%    Notes:
%     - 'bl' is assumed to be 'blue' as Matlab does
%
%    Examples:
%     invertcolor('rgb')
%     invertcolor('red')
%     invertcolor(['wh'; 'k'; 'bl'; 'y'])
%     invertcolor([0.5 0.5 0.5])
%     invertcolor([0.5 0.5 0.5],true)
%     invertcolor(jet)
%
%    See also: COLORMAP

%     Version History:
%        May   4, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   4, 2010 at 19:15 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% default flag
if(nargin==1 || isempty(flag)); flag=false; end

% check flag
if(~isscalar(flag) || (~isnumeric(flag) && ~islogical(flag)))
    error('seizmo:invertcolor:badFlag',...
        'FLAG must be TRUE or FALSE!');
end

% valid short names
sn='bcgkmrwy';
sni='yrmwgckb';

% cell string of long+short names
ln={'b' 'c' 'g' 'k' 'm' 'r' 'w' 'y' ...
    'blue' 'cyan' 'green' 'black' 'magenta' 'red' 'white' 'yellow'};
lni={'y' 'r' 'm' 'w' 'g' 'c' 'k' 'b' ...
    'yellow' 'red' 'magenta' 'white' 'green' 'cyan' 'black' 'blue'};

% allow Nx3 (must be between 0 & 1)
if(isreal(c) && size(c,2)==3 && all(c(:)<=1 & c(:)>=0))
    % color triplet
    if(flag)
        % high contrast mode
        c=1-round(c);
    else
        c=1-c;
    end
elseif(ischar(c))
    % if any chars outside short name list assume long name
    [tf,idx]=ismember(c(:),sn(:));
    if(all(tf))
        % short names
        c(:)=sni(idx);
    else
        % fall back to long/short names in rows
        % - using strmatch so we catch partial names
        % - in case of blue/black conflict, choose blue like matlab
        c=cellstr(c);
        for i=1:numel(c)
            idx=strmatch(c{i},ln);
            if(isempty(idx))
                error('seizmo:invertcolor:badColor',...
                    'Unknown Color: %s',c{:});
            end
            c(i)=lni(idx(1));
        end
    end
elseif(iscellstr(c))
    %long/short names in cells
    % - using strmatch so we catch partial names
    % - in case of blue/black conflict, choose blue like matlab
    for i=1:numel(c)
        idx=strmatch(c{i},ln);
        if(isempty(idx))
            error('seizmo:invertcolor:badColor',...
                'Unknown Color: %s',c{:});
        end
        c(i)=lni(idx(1));
    end
end

end
