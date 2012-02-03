function [idx]=vertcmp(data)
%VERTCMP    Finds vertical component SEIZMO records
%
%    Usage:    idx=vertcmp(data)
%
%    Description:
%     IDX=VERTCMP(DATA) returns the indices in matrix IDX for the records
%     in SEIZMO struct DATA that are vertically oriented (CMPINC is 0 or
%     180).  Special checks are made for inclined components (CMPINC is not
%     0, 90, or 180), vertically oriented North/East channels, horizontally
%     oriented Z channels, multiple vertical components for the same
%     stream, and varying azimuth or inclination for a component.
%
%    Notes:
%     - Any streams with reported problems will not be in the returned
%       results.  Think of it as encouragement to fix your data.
%
%    Examples:
%     % Vertical components in dataset without suspect orientations:
%     vdata=data(vertcmp(data));
%
%    See also: HORZCMP, HORZPAIRS, ROTATE

%     Version History:
%        Sep. 29, 2010 - initial version
%        Nov. 30, 2011 - doc update, fixed comments, all checks on all
%                        streams, remove 2nd output
%        Dec. 21, 2011 - forgot checking data structure and general header
%        Jan. 28, 2012 - pass char array to strnlen
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2012 at 16:00 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt to find verticals
try
    % check header
    data=checkheader(data);
    
    % get cmpinc, cmpaz, kname
    [cmpinc,cmpaz,kname]=getheader(data,'cmpinc','cmpaz','kname');
    
    % form stream name and component name
    sname=strcat(kname(:,1),'.',kname(:,2),'.',kname(:,3),'.',...
        strnlen(char(kname(:,4)),2));
    cname=kname(:,4);
    
    % extract the third letter of the kcmpnm
    cmpstr=char(cname); cmpstr=cellstr(cmpstr(:,3));

    % horizontals logical
    ishorz=cmpinc==90;

    % record indices
    idx=(1:numel(data))';

    % get stream index
    [streamcode,tmp,idx2]=unique(sname);

    % loop over streams
    isok=false(size(idx2));
    idx3=nan(size(idx2)); 
    for i=1:max(idx2)
        % start off as stream is ok
        bad=false;
        
        % first lets check out the non-horizontals for
        % potential metadata failures
        stream=(idx2==i & ~ishorz);
        
        % check for vertical N,E
        if(any(strcmpi(cmpstr(stream),'N')) ...
                || any(strcmpi(cmpstr(stream),'E')))
            warning('seizmo:vertcmp:verticalHorizontal',...
                ['Record(s):\n' sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name) ...
                '\nNorth/East components oriented vertically!']);
            bad=true;
        end
        
        % check for inclined anything
        if(any(cmpinc(stream)~=0 & cmpinc(stream)~=180))
            warning('seizmo:vertcmp:verticalHorizontal',...
                ['Record(s):\n' sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name) ...
                '\nComponent(s) oriented at strange incline!']);
            bad=true;
        end
        
        % check for multiple vertical components in stream
        if(numel(unique(cmpstr(stream)))>1)
            warning('seizmo:vertcmp:multipleVerticals',...
                ['Record(s):\n' sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name) ...
                '\nMultiple vertical components!']);
            bad=true;
        end
        
        % check that there is only one azimuth for the verticals
        if(numel(unique(cmpaz(stream)))>1)
            warning('seizmo:vertcmp:multiAziVerticals',...
                ['Record(s):\n' ...
                sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name)...
                '\nAzimuth of vertical component not consistent!']);
            bad=true;
        end
        
        % now on to the horizontals
        stream=(idx2==i & ishorz);
        
        % check for horizontal Z
        if(any(strcmpi(cmpstr(stream),'Z')))
            %tmp=stream(strcmpi(cmpstr(stream),'Z'));
            warning('seizmo:vertcmp:horizontalVertical',...
                ['Record(s):\n' sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name) ...
                '\nZ component oriented horizontally!']);
            bad=true;
        end
        
        % now on to either
        stream=(idx2==i);
        
        % get component index
        [componentcode,tmp,idx3(stream)]=unique(cname(stream));

        % get number of components
        ncmp=max(idx3(stream));

        % check that there is only one inclination per component
        for j=1:ncmp
            if(numel(unique(cmpinc(stream & idx3==j)))>1)
                warning('seizmo:vertcmp:multiOrientCmp',...
                    ['Record(s):\n' ...
                    sprintf('%d ',find(stream & idx3==j)) ...
                    '\nFilename(s):\n' ...
                    sprintf('%s\n',data(stream & idx3==j).name)...
                    '\nInclination for component not consistent!']);
                bad=true;
            end
        end

        % skip if anything found wrong above
        if(bad)
            continue;
        else
            % all good!
            isok(stream)=true;
        end
    end

    % return indices of verticals from good streams
    idx=idx(isok & ~ishorz);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
