function [rgb]=name2rgb(name)
%NAME2RGB    Converts short/long color names to RGB triplets
%
%    Usage:    rgb=name2rgb(name)
%
%    Description: RGB=NAME2RGB(NAME) converts short/long color names to
%     red-green-blue triplet values.  NAME must be a string or a cell array
%     of strings.  RGB is a Nx3 array of triplets corresponding to the
%     color names in NAME.
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
%     Make a colormap with all valid names:
%      colormap(name2rgb('roylgacsbvmpkw'));
%
%    See also: INVERTCOLOR, COLORMAP

%     Version History:
%        May  11, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  11, 2010 at 12:00 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% valid short names
sn='roylgacsbvmpkw';
sn_rgb=[1 0 0; 1 .5 0; 1 1 0; .5 1 0; 0 1 0; 0 1 .5; 0 1 1; 0 .5 1; 0 0 1;
    .5 0 1; 1 0 1; 1 0 .5; 0 0 0; 1 1 1];

% cell string of long+short names
ln={'r' 'o' 'y' 'l' 'g' 'a' 'c' 's' 'b' 'v' 'm' 'p' 'k' 'w' ...
    'red' 'orange' 'yellow' 'limegreen' 'green' 'aquamarine' 'cyan' ...
    'skyblue' 'blue' 'violet' 'magenta' 'pink' 'black' 'white'};
ln_rgb=[sn_rgb; sn_rgb];

% allow char/cellstr
if(ischar(name))
    % if any chars outside short name list assume long name
    [tf,idx]=ismember(name(:),sn(:));
    if(all(tf))
        % short names
        rgb=sn_rgb(idx,:);
    else
        % fall back to long/short names in rows
        % - using strmatch so we catch partial names
        % - in case of blue/black conflict, choose blue like matlab
        name=cellstr(name);
        rgb=nan(numel(name),3);
        for i=1:numel(name)
            idx=strmatch(name{i},ln);
            if(isempty(idx))
                error('seizmo:name2rgb:badName',...
                    'Unknown Color: %s',name{:});
            end
            rgb(i,:)=ln_rgb(idx(1),:);
        end
    end
elseif(iscellstr(name))
    %long/short names in cells
    % - using strmatch so we catch partial names
    % - in case of blue/black conflict, choose blue like matlab
    rgb=nan(numel(name),3);
    for i=1:numel(name)
        idx=strmatch(name{i},ln);
        if(isempty(idx))
            error('seizmo:name2rgb:badName',...
                'Unknown Color: %s',name{:});
        end
        rgb(i,:)=ln_rgb(idx(1),:);
    end
else
    error('seizmo:name2rgb:badName',...
        'NAME must be a string or a cell array of strings!');
end

end
