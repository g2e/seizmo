function [c]=invertcolor(c,flag)
%INVERTCOLOR    Inverts colors given as rgb triplet or as short/long names
%
%    Usage:    ci=invertcolor(c)
%              ci=invertcolor(c,hcmode)
%
%    Description:
%     CI=INVERTCOLOR(C) returns the "opposite" rgb color values CI to those
%     in C.  This is essentially CI=1-C but with special handling of short
%     and long color names (where it returns the opposing short or long
%     color name).  Note that "gray" [0.5 0.5 0.5] here is the opposite of
%     itself and so will not provide any contrast (see next usage form for
%     high contrast colors).
%
%     CI=INVERTCOLOR(C,HCMODE) toggles high contrast mode.  HCMODE set to
%     TRUE will find the color CI with the highest contrast to C using the
%     algorithm CI=1-ROUND(C).  This provides high contrast for grays
%     as all colors returned in high contrast mode are fully (de)saturated.
%     The default for FLAG is FALSE (return rgb color inverse).
%
%    Notes:
%     - 'bl' is assumed to be 'blue' as Matlab does
%     - Supports an extended set of color names:
%       SHORTNAME   LONGNAME         RGB
%       r           red              [1   0   0  ]
%       o           orange           [1   0.5 0  ]
%       y           yellow           [1   1   0  ]
%       l           limegreen        [0.5 1   0  ]
%       g           green            [0   1   0  ]
%       a           aquamarine       [0   1   0.5]
%       c           cyan             [0   1   1  ]
%       s           skyblue          [0   0.5 1  ]
%       b           blue             [0   0   1  ]
%       v           violet           [0.5 0   1  ]
%       m           magenta          [1   0   1  ]
%       p           pink             [1   0   0.5]
%       k           black            [0   0   0  ]
%       w           white            [1   1   1  ]
%
%    Examples:
%     invertcolor('rgb')
%     invertcolor('red')
%     invertcolor(['wh'; 'k'; 'bl'; 'y'])
%     invertcolor([0.5 0.5 0.5])
%     invertcolor([0.5 0.5 0.5],true)
%     invertcolor(jet)
%
%    See also: COLORMAP, NAME2RGB

%     Version History:
%        May   4, 2010 - initial version
%        May  11, 2010 - added support for the extended color name set
%        Aug.  4, 2010 - see also NAME2RGB
%        Sep. 13, 2010 - remove see also itself
%        Feb.  1, 2011 - renamed flag to hcmode for clarity
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 14:15 GMT

% todo:
% - hsv support

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
sn='roylgacsbvmpkw';
sni='csbvmproylgawk';

% cell string of long+short names
ln={'r' 'o' 'y' 'l' 'g' 'a' 'c' 's' 'b' 'v' 'm' 'p' 'k' 'w' ...
    'red' 'orange' 'yellow' 'limegreen' 'green' 'aquamarine' 'cyan' ...
    'skyblue' 'blue' 'violet' 'magenta' 'pink' 'black' 'white'};
lni={'c' 's' 'b' 'v' 'm' 'p' 'r' 'o' 'y' 'l' 'g' 'a' 'w' 'k' ...
    'cyan' 'skyblue' 'blue' 'violet' 'magenta' 'pink' 'red' 'orange' ...
    'yellow' 'limegreen' 'green' 'aquamarine' 'white' 'black'};

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
