function [idx1,idx2]=vertcmp(data)
%VERTCMP    Finds vertical component SEIZMO records
%
%    Usage:    [idx1,idx2]=vertcmp(data)
%
%    Description:
%     [IDX1,IDX2]=VERTCMP(DATA) returns the indices in matrix IDX1 for the
%     records in SEIZMO struct DATA that are vertically oriented (CMPINC is
%     0 or 180).  Special checks are made for inclined components (CMPINC
%     is not 0, 90, or 180), vertically oriented North/East channels,
%     multiple vertical components for the same stream, and varying azimuth
%     or inclination.  IDX2 contains the corresponding stream index.
%
%    Notes:
%
%    Examples:
%     % Vertical components in dataset without suspect orientations:
%     vdata=data(vertcmp(data));
%
%    See also: HORZPAIRS, ROTATE

%     Version History:
%        Sep. 29, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 29, 2010 at 16:00 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt pairing of horizontals
try
    % get cmpinc, cmpaz, kname
    [cmpinc,cmpaz,kname]=getheader(data,'cmpinc','cmpaz','kname');
    
    % form stream name and component name
    sname=strcat(kname(:,1),'.',kname(:,2),'.',kname(:,3),'.',...
        strnlen(kname(:,4),2));
    cname=kname(:,4);
    
    % extract the third letter of the kcmpnm
    cmpstr=char(cname); cmpstr=cellstr(cmpstr(:,3));

    % only horizontals
    ishorz=cmpinc==90;

    % record indices
    idx1=(1:numel(data))';

    % get stream index
    [streamcode,tmp,idx2]=unique(sname);

    % loop over streams
    offset=0; % offset to assure pair indices do not skip numbers
    isok=false(size(idx2));
    idx3=nan(size(idx2)); 
    for i=1:max(idx2)
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
            offset=offset+1;
            continue;
        end
        
        % check for inclined anything
        if(any(cmpinc(stream)~=0 & cmpinc(stream)~=180))
            warning('seizmo:vertcmp:verticalHorizontal',...
                ['Record(s):\n' sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name) ...
                '\nComponent(s) oriented at strange incline!']);
            offset=offset+1;
            continue;
        end
        
        % check for multiple vertical components in stream
        if(numel(unique(cmpstr(stream)))>1)
            warning('seizmo:vertcmp:multipleVerticals',...
                ['Record(s):\n' sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name) ...
                '\nMultiple vertical components!']);
            offset=offset+1;
            continue;
        end
        
        % check that there is only one orientation per component
        if(numel(unique(cmpaz(stream)))>1)
            warning('seizmo:vertcmp:multiOrientCmp',...
                ['Record(s):\n' ...
                sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name)...
                '\nOrientation for component not consistent!']);
            offset=offset+1;
            continue;
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
            offset=offset+1;
            continue;
        end
        
        % now on to either
        stream=(idx2==i);
        
        % get component index
        [componentcode,tmp,idx3(stream)]=unique(cname(stream));

        % get number of horizontal components
        ncmp=max(idx3(stream));

        % check that there is only one inclination per component
        badcmp=false;
        for j=1:ncmp
            if(numel(unique(cmpinc(stream & idx3==j)))>1)
                warning('seizmo:vertcmp:multiOrientCmp',...
                    ['Record(s):\n' ...
                    sprintf('%d ',find(stream & idx3==j)) ...
                    '\nFilename(s):\n' ...
                    sprintf('%s\n',data(stream & idx3==j).name)...
                    '\nOrientation for single component not consistent!']);
                badcmp=true;
            end
        end

        % skip if cmp w/ 2+ inclinations
        if(badcmp); offset=offset+1; continue; end
        
        % update pair number
        idx2(stream)=i-offset;
        
        % is ok
        isok(stream)=true;
    end

    % return indices of horizontal pairs
    idx1=idx1(isok & ~ishorz);
    idx2=idx2(isok & ~ishorz);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
