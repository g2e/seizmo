function [ref1,ref2,offset1,offset2,fill,filler,trim]=cutparam(varargin)
%CUTPARAM    Parses inputs defining the cut window for CUTIM/RPDW
%
%    Description: CUTPARAM(LIST,OF,WINDOW,PARAMETERS) parses the parameters
%     passed to CUTIM and RPDW passing the results back.  This includes the
%     those parameters that define the window as well as other options 
%     (fill etc).  Parameters go through a number of basic checks (but by
%     no means complete).
%
%    Notes:
%     - supplying empty arguments can have unexpected results
%
%    System requirements: Matlab 7
%
%    Data requirements: Parameters should be either strings or numeric.
%
%    Header changes: NONE
%
%    Usage: [ref1,ref2,offset1,offset2,fill,filler,trim]=cutparam(args)
%
%    See also: cutim, rpdw

%     Version History:
%        Apr. 17, 2008 - initial version
%        Apr. 18, 2008 - bugfix
%        June 24, 2008 - added more checks and a major documentation update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 24, 2008 at 00:15 GMT

% defaults
trim=true;
fill=false;
filler=0;
ref1='b';
ref2='e';
offset1=0;
offset2=0;

% read in cut parameters
skip=false;     % skip next argument (already read in)
ref1unset=true; % is REF1/OFFSET1 not given
ref2unset=true; % is REF2/OFFSET2 not given
for i=1:nargin
    % skip flag
    if(skip); skip=false; continue; end
    
    % skip empty
    if(isempty(varargin{i})); continue; end
    
    % offset without reference
    if(isnumeric(varargin{i}))
        % first offset (set reference to 'z')
        if(ref1unset)
            ref1='z';
            offset1=varargin{i};
            ref1unset=false;
        % second offset (use first reference)
        elseif(ref2unset)
            ref2=ref1;
            offset2=varargin{i};
            ref2unset=false;
        % only two offsets allowed
        else
            error('SAClab:cutparam:badInput','Too many offsets!')
        end
    % string (could be reference or other option)
    elseif(ischar(varargin{i}))
        skip=true;
        
        % make sure value is scalar if an option
        if(any(strcmpi({'fill' 'filler' 'trim'},varargin{i})))
            if(isempty(varargin{i+1}))
                % leave as default by skipping
                continue;
            elseif(~isscalar(varargin{i+1}) || ...
                    (~isnumeric(varargin{i+1}) && ...
                    ~islogical(varargin{i+1})))
                error('SAClab:cutparam:badInput',...
                    'Option value must be a numeric or logical scalar')
            end
        end
        
        % option fill
        if(strcmpi('fill',varargin{i}))
            fill=varargin{i+1};
        % option filler
        elseif(strcmpi('filler',varargin{i}))
            if(islogical(varargin{i+1}))
                error('SAClab:cutparam:badInput',...
                    'Filler value must be a numeric scalar');
            end
            filler=varargin{i+1};
        % option trim
        elseif(strcmpi('trim',varargin{i}))
            trim=varargin{i+1};
        % otherwise assume reference
        else
            % first reference
            if(ref1unset)
                ref1=varargin{i};
                ref1unset=false;
                
                % first offset
                if(isempty(varargin{i+1}))
                    continue;
                elseif(isnumeric(varargin{i+1}))
                    offset1=varargin{i+1};
                else
                    offset1=0;
                    skip=false;
                end
            % second reference
            elseif(ref2unset)
                ref2=varargin{i};
                ref2unset=false;
                
                % second offset
                if(isempty(varargin{i+1}))
                    continue;
                elseif(isnumeric(varargin{i+1}))
                    offset2=varargin{i+1};
                else
                    offset2=0;
                    skip=false;
                end
            % only two references allowed
            else
                error('SAClab:cutparam:badInput','Too many references!')
            end
        end
    else
        error('SAClab:cutparam:badInput',...
            'Window parameters must be strings or numeric!')
    end
end

end
