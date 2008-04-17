function [ref1,ref2,offset1,offset2,fill,filler,trim]=cutparam(varargin)
%CUTPARAM    Parses inputs defining the cut window

% defaults
trim=true;
fill=false;
filler=0;
ref1='b';
ref2='e';
offset1=0;
offset2=0;

% read in cut parameters
skip=false;
ref1unset=true;
ref2unset=true;
for i=1:nargin
    % skip flag
    if(skip); skip=false; continue; end
    
    % offset without reference
    if(isnumeric(varargin{i}))
        % first offset (set reference to 'z')
        if(ref1unset)
            ref1='z';
            offset1=varargin{i};
        % second offset (use first reference)
        elseif(ref2unset)
            ref2=ref1;
            offset2=varargin{i};
        % only two offsets allowed
        else
            error('SAClab:cutim:badInput','Too many offsets')
        end
    % string (could be reference or other option)
    else
        skip=true;
        % option fill
        if(strcmpi('fill',varargin{i}))
            fill=varargin{i+1};
        % option filler
        elseif(strcmpi('filler',varargin{i}))
            filler=varargin{i+1};
        % option removedataless
        elseif(strcmpi('removedataless',varargin{i}))
            trim=varargin{i+1};
        % otherwise assume reference
        else
            % first reference
            if(ref1unset)
                ref1=varargin{i};
                ref1unset=false;
                
                % first offset
                if(isnumeric(varargin{i+1}))
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
                if(isnumeric(varargin{i+1}))
                    offset2=varargin{i+1};
                else
                    offset2=0;
                    skip=false;
                end
            % only two references allowed
            else
                error('SAClab:cutim:badInput','Too many references!')
            end
        end
    end
end

end

