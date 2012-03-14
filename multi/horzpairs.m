function [idx1,idx2,idx3]=horzpairs(data,varargin)
%HORZPAIRS    Returns indice arrays for pairing horizontal SEIZMO records
%
%    Usage:    [idx1,idx2,idx3]=horzpairs(data)
%              horzpairs(...,'requiredcharfields',fields,...)
%              horzpairs(...,'requiredrealfields',fields,...)
%
%    Description:
%     [IDX1,IDX2,IDX3]=HORZPAIRS(DATA) pairs orthogonal horizontal
%     components for records in SEIZMO struct DATA.  Records are paired
%     based on their KNETWK, KSTNM, KHOLE, KCMPNM, LEVEN, NCMP, CMPINC, &
%     CMPAZ header fields.  Pairs are identified first by having the same
%     stream (see GETSTREAMIDX) + LEVEN + NCMP + CMPINC==90 and then by
%     azimuth.  Individual components must have 1 orientation only to be
%     paired.  Streams with more than 2 different horizontals are not
%     paired.  Warnings are issued for such cases as well as when a
%     horizontal pair is non-orthogonal, a horizontal (N,E) is oriented
%     non-horizontally, a vertical (Z) is oriented horizontally and when
%     any component is not vertical or horizontal.  Anytime a warning is
%     issued for a stream, that stream will not be paired (fix your
%     metadata folks).  IDX1 gives the indices of records in DATA that have
%     been paired.  IDX2 gives the corresponding pair indices (run
%     max(IDX2) to get the number of pairs).  IDX3 gives the component
%     indices (always 1 or 2).
%
%     [IDX1,IDX2,IDX3]=HORZPAIRS(...,'REQUIREDCHARFIELDS',FIELDS,...)
%     allows changing the character fields required to be equal between
%     records before they are paired.  The list is a cellstring array.  The
%     default is: {}.  Note that pairing based on KNETWK, KSTNM, KHOLE,
%     KCMPNM header fields is hard-coded.
%
%     [IDX1,IDX2,IDX3]=HORZPAIRS(...,'REQUIREDREALFIELDS',FIELDS,...)
%     changes the numerical fields required to be equal between records
%     before they are paired.  The list must be a cellstring array.  The
%     default is: {'delta'}.  Note that CMPINC, LEVEN and NCMP are also
%     required but cannot be removed from the list.
%
%    Notes:
%     - Run FIXDELTA first to take care of small differences in sample
%       rates caused by floating point inaccuracies!
%
%    Examples:
%     % Get pairs and loop over pairs:
%     [idx1,idx2,idx3]=horzpairs(data);
%     for i=1:max(idx2)
%         % record indices for this pair
%         ridx=idx1(idx2==i);
%
%         % separate indices based on component
%         cidx1=idx1(idx2==i & idx3==1);
%         cidx2=idx1(idx2==i & idx3==2);
%
%         ... % do something useful here % ...
%     end
%
%    See also: ROTATE, ISORTHOGONAL, GETSTREAMIDX, GETCOMPONENTIDX

%     Version History:
%        Feb. 22, 2010 - initial version
%        Feb. 24, 2010 - more checks, 1st cmp lags 2nd by 90deg
%        Sep. 29, 2010 - add filename output to warnings, more warnings for
%                        verticals, check for multiple inclinations
%        Feb. 11, 2011 - dropped versioninfo caching
%        Jan. 28, 2012 - doc update, char to strnlen, drop SEIZMO global
%        Feb.  7, 2012 - fixed a couple warning ids
%        Mar. 13, 2012 - use getheader improvements
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2012 at 15:05 GMT

% todo:

% check nargin
if(mod(nargin-1,2))
    error('seizmo:horzpairs:badNumInputs',...
        'Bad number of arguments!');
end

% check headers
data=checkheader(data);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt pairing of horizontals
try
    % defaults
    option.REQUIREDCHARFIELDS={};
    option.REQUIREDREALFIELDS={'delta'};

    % get options from command line
    for i=1:2:nargin-1
        if(~ischar(varargin{i}))
            error('seizmo:horzpairs:badInput',...
                'Options must be specified as a string!');
        end
        if(~isempty(varargin{i+1}))
            option.(upper(varargin{i}))=varargin{i+1};
        end
    end
    
    % check options
    fields=fieldnames(option);
    for i=1:numel(fields)
        % get value of field and do a basic check
        value=option.(fields{i});

        % specific checks
        switch lower(fields{i})
            case {'requiredcharfields' 'requiredrealfields'}
                % fix char arrays
                if(ischar(value))
                    value=cellstr(value);
                    option.(fields{i})=value;
                end
                if(~iscellstr(value))
                    error('seizmo:horzpairs:badInput',...
                        '%s option must be a cellstr of header fields!',...
                        fields{i});
                end
            otherwise
                error('seizmo:horzpairs:badInput',...
                    'Unknown option: %s !',fields{i});
        end
    end
    
    % get cmpinc, cmpaz, kname
    [cmpinc,cmpaz,kname,ncmp,leven]=getheader(data,...
        'cmpinc','cmpaz','kname','ncmp','leven lgc');
    
    % get soft requirements
    szreal=size(option.REQUIREDREALFIELDS); reqreal=cell(szreal);
    szchar=size(option.REQUIREDCHARFIELDS); reqchar=cell(szchar);
    if(prod(szreal)~=0)
        [reqreal{:}]=getheader(data,option.REQUIREDREALFIELDS{:});
    end
    if(prod(szchar)~=0)
        [reqchar{:}]=getheader(data,option.REQUIREDCHARFIELDS{:});
    end
    
    % get stream+req name and component name
    sname=strcat(kname(:,1),'.',kname(:,2),'.',kname(:,3),'.',...
        strnlen(char(kname(:,4)),2),'_',strcat('',reqchar{:}),'_',...
        strcat('',reqreal{:}),'_',leven,'_',num2str(ncmp));
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
    idx3=nan(size(idx2)); ispair=false(size(idx2));
    for i=1:max(idx2)
        % first lets check out the non-horizontals for
        % potential metadata failures
        stream=(idx2==i & ~ishorz);
        
        % check for vertical N,E
        if(any(strcmpi(cmpstr(stream),'N')) ...
                || any(strcmpi(cmpstr(stream),'E')))
            %tmp=stream(strcmpi(cmpstr(stream),'N') ...
            %    | strcmpi(cmpstr(stream),'E'));
            warning('seizmo:horzpairs:verticalHorizontal',...
                ['Record(s):\n' sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name) ...
                '\nNorth/East components oriented vertically!']);
            offset=offset+1;
            continue;
        end
        
        % check for inclined anything
        if(any(cmpinc(stream)~=0 & cmpinc(stream)~=180))
            %tmp=stream(cmpinc(stream)~=0 & cmpinc(stream)~=180);
            warning('seizmo:horzpairs:verticalHorizontal',...
                ['Record(s):\n' sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name) ...
                '\nComponent(s) oriented at strange incline!']);
            offset=offset+1;
            continue;
        end
        
        % check for multiple vertical components in stream
        if(numel(unique(cmpstr(stream)))>1)
            warning('seizmo:horzpairs:multipleVerticals',...
                ['Record(s):\n' sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name) ...
                '\nMultiple vertical components!']);
            offset=offset+1;
            continue;
        end
        
        % check that there is only one orientation per component
        if(numel(unique(cmpaz(stream)))>1)
            warning('seizmo:horzpairs:multiOrientCmp',...
                ['Record(s):\n' ...
                sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name)...
                '\nOrientation for component not consistent!']);
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
                warning('seizmo:horzpairs:multiOrientCmp',...
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
        
        % now on to the horizontals
        stream=(idx2==i & ishorz);
        
        % skip if none
        if(sum(stream)==0)
            offset=offset+1;
            continue;
        end
        
        % check for horizontal Z
        if(any(strcmpi(cmpstr(stream),'Z')))
            %tmp=stream(strcmpi(cmpstr(stream),'Z'));
            warning('seizmo:horzpairs:horizontalVertical',...
                ['Record(s):\n' sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name) ...
                '\nZ component oriented horizontally!']);
            offset=offset+1;
            continue;
        end

        % get component index
        [componentcode,tmp,idx3(stream)]=unique(cname(stream));

        % get number of horizontal components
        ncmp=max(idx3(stream));

        % check that there is only one orientation per component
        badcmp=false;
        for j=1:ncmp
            if(numel(unique(cmpaz(stream & idx3==j)))>1)
                warning('seizmo:horzpairs:multiOrientCmp',...
                    ['Record(s):\n' ...
                    sprintf('%d ',find(stream & idx3==j)) ...
                    '\nFilename(s):\n' ...
                    sprintf('%s\n',data(stream & idx3==j).name)...
                    '\nOrientation for single component not consistent!']);
                badcmp=true;
            end
        end

        % skip if only one horizontal or cmp w/ 2+ orientations
        if(ncmp==1 || badcmp); offset=offset+1; continue; end

        % error if more than 2 horizontals
        if(ncmp>2)
            warning('seizmo:horzpairs:multiHorzCmp',...
                ['Record(s):\n' sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name) ...
                '\nMore than 2 horizontal components for this station!']);
            offset=offset+1;
            continue;
        end
        
        % error if only one orientation
        az=unique(cmpaz(stream));
        if(numel(az)==1)
            warning('seizmo:horzpairs:nonOrthoHorz',...
                ['Record(s):\n' sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name) ...
                '\nAzimuth(s): ' sprintf('%g ',az) ...
                '\nHorizontal pair shares the same azimuth!']);
            offset=offset+1;
            continue;
        end

        % now check that pair is orthogonal
        if(~isorthogonal([90 az(1)],[90 az(2)]))
            warning('seizmo:horzpairs:nonOrthoHorz',...
                ['Record(s):\n' sprintf('%d ',find(stream)) ...
                '\nFilename(s):\n' sprintf('%s\n',data(stream).name) ...
                '\nAzimuth(s): ' sprintf('%g ',az) ...
                '\nHorizontal pair is not orthogonal!']);
            offset=offset+1;
            continue;
        else
            ispair(stream)=true;
        end
        
        % update pair number
        idx2(stream)=i-offset;
        
        % force second cmp to lead by 90deg
        % - this keeps things predictable
        cmp1=stream & idx3==1;
        cmp2=stream & idx3==2;
        az=[cmpaz(find(cmp1,1)) cmpaz(find(cmp2,1))];
        if(abs(mod(az(1)-90,360)-az(2))<1)
            % switch cmp indices
            idx3(cmp1)=2;
            idx3(cmp2)=1;
        end
    end

    % return indices of horizontal pairs
    idx1=idx1(ispair & ishorz);
    idx2=idx2(ispair & ishorz);
    idx3=idx3(ispair & ishorz);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
