function [gmt]=readgmt(file,type,marker,ll)
%READGMT    Reads the GMT vector file-format (ascii)
%
%    Usage:    gmt=readgmt(file)
%              gmt=readgmt(file,type)
%              gmt=readgmt(file,type,marker)
%
%    Description:
%     GMT=READGMT(FILE) reads in a GMT ascii data file (sometimes called
%     the GMT vector ascii file).  These files are almost always denoted
%     with a .gmt extension.  The contents (usually 1 or more lines) is
%     output as the struct GMT which contains the following fields:
%      .type        - segment type (default is 'line')
%      .comments    - full line comments before/in this segment
%      .header      - text on segment's header line (denoted with a '>')
%      .latitude    - latitude position(s) of this segment
%      .longitude   - longitude position(s) of this segment
%      .text        - remaining text on lon/lat line(s) of this segment
%
%     GMT=READGMT(FILE,TYPE) alters the feature type to TYPE.  The default
%     is 'line' and is currently unused.
%
%     GMT=READGMT(FILE,TYPE,MARKER) alters the segment separator (aka
%     marker) to MARKER, a single ascii character.  The default is '>'.
%
%     GMT=READGMT(FILE,TYPE,MARKER,LATLON) sets if the points are given as
%     [lat lon] or [lon lat].  The default is FALSE ([lon lat]).
%
%    Notes:
%     - Currently only handles the old GMT format.  OGM extensions are not
%       handled yet (ie they are ignored).
%
%    Examples:
%     % Read plate boundary .gmt file:
%     file=[fileparts(which('readgmt.m')) filesep 'plate_boundaries.gmt'];
%     plates=readgmt(file);
%
%    See also: WRITEGMT

%     Version History:
%        Jan. 20, 2011 - initial version
%        Jan. 24, 2011 - support comma delimited
%        Jan. 31, 2011 - latlon arg
%        Jan. 27, 2014 - abs path exist fix
%        Feb.  9, 2014 - use readtxt
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2014 at 10:35 GMT

% todo

% check nargin
error(nargchk(0,3,nargin));

% default type/marker
if(nargin<2 || isempty(type)); type=[]; end
if(nargin<3 || isempty(marker)); marker='>'; end
if(nargin<4 || isempty(ll)); ll=false; end

% check type/marker
validtype={...
    'point' 'multipoint' ...
    'line' 'multiline' ...
    'polygon' 'multipolygon'};
if(~isempty(type) && (~isstring(type) || ~any(strcmpi(type,validtype))))
    error('seizmo:readgmt:badInput',...
        ['TYPE must be one of the following:\n' ...
        sprintf('''%s'' ',validtype{:})]);
elseif(~ischar(marker) || ~isscalar(marker))
    error('seizmo:readgmt:badInput',...
        'MARKER must be a single character (like ''>'')!');
elseif(~isscalar(ll) || (~islogical(ll) && ~isreal(ll)))
    error('seizmo:readgmt:badInput',...
        'LATLON must be either TRUE or FALSE!');
end

% read in file
filterspec={
    '*.gmt;*.GMT' 'GMT Files (*.gmt,*.GMT)';
    '*.txt;*.TXT' 'TXT Files (*.txt,*.TXT)';
    '*.*' 'All Files (*.*)'};
lines=getwords(readtxt(file,filterspec),sprintf('\n'));

% initialize output
gmt([])=struct('type',[],'comments',[],'header',[],...
    'latitude',[],'longitude',[],'text',[]);

% keep processing until through all lines
a=1; obj=1; nc=0; d=0;
nlines=numel(lines);
while(a<=nlines)
    % skip line if blank
    if(isempty(lines{a}))
        a=a+1;
        continue;
    end
    
    % process line
    words=getwords(lines{a});
    
    % skip line if blank
    if(isempty(words))
        a=a+1;
        continue;
    end
    
    % record comments
    % - we need to handle multi-segment comments
    %   - requires moving comments to next segment
    if(strcmp(words{1}(1),'#'))
        % remove comment character
        words{1}=words{1}(2:end);
        if(isempty(words{1})); words(1)=[]; end
        nc=nc+1;
        gmt(obj).comment{nc,1}=joinwords(words);
        a=a+1;
        continue;
    end
    
    % if marker encountered => new line/polygon
    if(strcmp(words{1}(1),marker))
        % set type
        if(isempty(type)); type='line'; end
        if(~strcmp(type,'line'))
            error('seizmo:readgmt:badGMT',...
                ['GMT File: %s\nLine %d: %s\n'...
                'File must not have multiple data types!'],...
                file,a,lines{a});
        end
        gmt(obj).type=type;
        
        % extract header
        words{1}=words{1}(2:end);
        if(isempty(words{1})); words(1)=[]; end
        gmt(obj).header=joinwords(words);
        
        % read in data until next segment or end
        a=a+1; d=0; nextsegment=false; csd=0;
        while(a<=nlines && ~nextsegment)
            % skip line if blank
            if(isempty(lines{a}))
                a=a+1;
                continue;
            end
            
            % process line
            words=getwords(lines{a});
            
            % skip line if blank
            if(isempty(words))
                a=a+1;
                continue;
            end
            
            % record comments
            if(strcmp(words{1}(1),'#'))
                % remove comment character
                words{1}=words{1}(2:end);
                if(isempty(words{1})); words(1)=[]; end
                nc=nc+1;
                csd=csd+1;
                gmt(obj).comment{nc,1}=joinwords(words);
                a=a+1;
                continue;
            end
            
            % data or new segment?
            if(strcmp(words{1}(1),marker))
                nextsegment=true;
            else
                % attempt splitting 1st word by commas
                words=[getwords(words{1},',') words(2:end)];
                
                % check has at least 2 numbers
                % - does not work with complex dd:mm:ss like positions
                if(numel(words)<2)
                    error('seizmo:readgmt:badData',...
                        ['GMT File: %s\nLine %d: %s\n'...
                        'Line must have at least LON LAT!'],...
                        file,a,lines{a});
                elseif(isnan(str2double(words{1})) ...
                        || isnan(str2double(words{1})))
                    error('seizmo:readgmt:badData',...
                        ['GMT File: %s\nLine %d: %s\n'...
                        'LON LAT must not be in DMS or have E/W/N/S!'],...
                        file,a,lines{a});
                end
                d=d+1;
                if(ll)
                    gmt(obj).latitude(d,1)=str2double(words{1});
                    gmt(obj).longitude(d,1)=str2double(words{2});
                else
                    gmt(obj).longitude(d,1)=str2double(words{1});
                    gmt(obj).latitude(d,1)=str2double(words{2});
                end
                gmt(obj).text{d,1}=joinwords(words(3:end));
                a=a+1;
                csd=0;
            end
        end
        
        % clean up segment (incrementing, moving comments)
        if(csd)
            gmt(obj+1).comment(1:csd,1)=gmt(obj).comment(end-csd+1:end,1);
            gmt(obj).comment(end-csd+1:end,1)=[];
            nc=csd;
        else
            nc=0;
        end
        obj=obj+1;
    else % hmmm...points then?
        % set type
        if(isempty(type)); type='point'; end
        if(~strcmp(type,'point'))
            error('seizmo:readgmt:badGMT',...
                ['GMT File: %s\nLine %d: %s\n'...
                'File must not have multiple data types!'],...
                file,a,lines{a});
        end
        gmt(obj).type=type;
        
        % attempt splitting 1st word by commas
        words=[getwords(words{1},',') words(2:end)];
        
        % check has at least 2 numbers
        % - does not work with complex dd:mm:ss like positions
        if(numel(words)<2)
            error('seizmo:readgmt:badData',...
                ['GMT File: %s\nLine %d: %s\n'...
                'Line must have at least LON LAT!'],...
                file,a,lines{a});
        elseif(isnan(str2double(words{1})) ...
                || isnan(str2double(words{1})))
            error('seizmo:readgmt:badData',...
                ['GMT File: %s\nLine %d: %s\n'...
                'LON LAT must not be in DMS or have E/W/N/S!'],...
                file,a,lines{a});
        end
        d=d+1;
        if(ll) % lat lon
            gmt(obj).latitude(d,1)=str2double(words{1});
            gmt(obj).longitude(d,1)=str2double(words{2});
        else % lon lat
            gmt(obj).longitude(d,1)=str2double(words{1});
            gmt(obj).latitude(d,1)=str2double(words{2});
        end
        gmt(obj).text{d,1}=joinwords(words(3:end));
        a=a+1;
        obj=obj+1;
    end
end

% delete last object if just a terminator
if(isempty(gmt(end).header) && isempty(gmt(end).latitude))
    gmt(end)=[];
end

end
