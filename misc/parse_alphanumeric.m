function [an,isnum]=parse_alphanumeric(str)
%PARSE_ALPHANUMERIC    Split alphanumeric string into words & numbers
%
%    Usage:    [an,isnum]=parse_alphanumeric(str)
%
%    Description:
%     [AN,ISNUM]=PARSE_ALPHANUMERIC(STR) parses out alphabet and digit
%     sequences from character string STR as a cell array of "words" and
%     numbers (converted to double) in AN.  ISNUM indicates the elements in
%     AN that are numeric.
%
%    Notes:
%
%    Examples:
%     % Split a date string:
%     parse_alphanumeric('2000may03')
%
%    See also: ISSTRPROP, GETWORDS, READTXT, VERSION_COMPARE

%     Version History:
%        Feb. 13, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 13, 2012 at 13:00 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% require character input
if(~ischar(str))
    error('seizmo:parse_alphanumeric:badInput',...
        'STR must be a character string!');
end

% identify alphanumeric bits
num=isstrprop(str,'digit');
alf=isstrprop(str,'alpha');

% loop over bits getting words & numbers
cnt=1; in=false; idx=nan(0,2); isnum=false(0,1);
for i=1:numel(str)
    if(num(i))
        % number
        if(in(cnt) && isnum(cnt))
            % continuing number
            continue;
        elseif(in(cnt) && ~isnum(cnt))
            % terminate alpha, begin number
            idx(cnt,2)=i-1;
            cnt=cnt+1;
            idx(cnt,1)=i;
            isnum(cnt)=true;
            in(cnt)=true;
        else
            % new number
            idx(cnt,1)=i;
            isnum(cnt)=true;
            in(cnt)=true;
        end
    elseif(alf(i))
        % word
        if(in(cnt) && isnum(cnt))
            % terminate number, begin word
            idx(cnt,2)=i-1;
            cnt=cnt+1;
            idx(cnt,1)=i;
            isnum(cnt)=false;
            in(cnt)=true;
        elseif(in(cnt) && ~isnum(cnt))
            % continuing word
            continue;
        else
            % new word
            idx(cnt,1)=i;
            isnum(cnt)=false;
            in(cnt)=true;
        end
    else
        % not either
        if(in(cnt) && isnum(cnt))
            % terminate number
            idx(cnt,2)=i-1;
            cnt=cnt+1;
            in(cnt)=false;
        elseif(in(cnt) && ~isnum(cnt))
            % terminate word
            idx(cnt,2)=i-1;
            cnt=cnt+1;
            in(cnt)=false;
        end
    end
end

% close final word/number
if(in(cnt)); idx(cnt,2)=i; end

% make cellstr
an=cell(1,numel(isnum));
for i=1:numel(isnum)
    if(isnum(i))
        an{i}=str2double(str(idx(i,1):idx(i,2)));
    else
        an{i}=str(idx(i,1):idx(i,2));
    end
end

end
