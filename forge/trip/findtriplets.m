function [idx1,idx2,idx3]=findtriplets(data,varargin)
%FINDTRIPLETS    Returns indexing for 3 component sets in a SEIZMO dataset
%
%    Usage:    [idx1,idx2,idx3]=findtriplets(data)
%              [...]=findtriplets(...,'requiredcharfields',fields,...)
%              [...]=findtriplets(...,'requiredrealfields',fields,...)
%
%    Description:
%     [IDX1,IDX2,IDX3]=FINDTRIPLETS(DATA) forms triplets of orthogonal sets
%     of components in SEIZMO struct DATA.  These triplets can then be used
%     to look at the data in 3D.  Triplets are found based on their KNETWK,
%     KSTNM, KHOLE, KCMPNM, LEVEN, NCMP, CMPINC, and CMPAZ header fields.
%     Additional header fields required to be equal can be set by passing
%     parameters described below with the default requiring the DELTA field
%     to be equal.  Warnings are issued for a variety of typical metadata
%     problems encountered for each potential triplet and only triplets
%     without metadata issues are returned (fix your metadata folks).  IDX1
%     gives the indices of records in DATA that are in a triplet.  IDX2
%     gives the corresponding triplet indices (run max(IDX2) to get the
%     number of triplets).  IDX3 gives the component indices (1, 2 or 3).
%
%     [...]=FINDTRIPLETS(...,'REQUIREDCHARFIELDS',FIELDS,...) changes the
%     character header fields required to be equal between triplets.  The
%     list is a cellstring array.  The default is: {}.  Note KNETWK, KSTNM,
%     KHOLE, KCMPNM header fields requirements are not optional.
%
%     [...]=FINDTRIPLETS(...,'REQUIREDREALFIELDS',FIELDS,...) changes the
%     numerical header fields required to be equal between triplets.  The
%     list must be a cellstring array.  The default is: {'delta'}.  Note
%     that LEVEN and NCMP are always required to be equal.
%
%    Notes:
%     - FINDTRIPLETS skips a triplet after the first metadata issue.  You
%       may have to run FINDTRIPLETS multiple times when fixing metadata
%       issues to resolve them all.
%
%    Examples:
%     % Use with no outputs to do metadata checking of your data:
%     findtriplets(data);
%
%    See also: MAKETRIPLETS, ROTATE3, PLOTPM3, PMMOVIE3, UNMAKETRIPLETS,
%              HORZPAIRS, HORZCMP, VERTCMP, ROTATE, PLOTPM, PMMOVIE

%     Version History:
%        July 30, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 30, 2012 at 15:05 GMT

% todo:

% check nargin
if(mod(nargin-1,2))
    error('seizmo:findtriplets:badNumInputs',...
        'Bad number of arguments!');
end

% check headers
data=checkheader(data);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt finding triplets
try
    % defaults
    option.REQUIREDCHARFIELDS={};
    option.REQUIREDREALFIELDS={'delta'};

    % get options from command line
    for i=1:2:nargin-1
        if(~ischar(varargin{i}))
            error('seizmo:findtriplets:badInput',...
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
                    error('seizmo:findtriplets:badInput',...
                        '%s option must be a cellstr of header fields!',...
                        fields{i});
                end
            otherwise
                error('seizmo:findtriplets:badInput',...
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

    % record indices
    idx1=(1:numel(data))';

    % get stream index
    [streamcode,tmp,idx2]=unique(sname);
    
    % loop over streams
    offset=0; % offset to assure trip indices do not skip numbers
    idx3=nan(size(idx2)); istriplet=false(size(idx2));
    for i=1:max(idx2)
        % records in stream
        stream=(idx2==i);
        
        % get component index
        [componentcode,idx4,idx3(stream)]=unique(cname(stream));
        
        % get number of components
        nstreamcmp=max(idx3(stream));
        
        % check each cmp is consistent
        badcmp=false;
        for j=1:nstreamcmp
            if(numel(unique(cmpaz(stream & idx3==j)))>1)
                warning('seizmo:findtriplets:cmpChangesAzimuth',...
                    ['Azimuth for component not consistent!\n' ...
                    'Filename(s):\n' ...
                    sprintf(' %s\n',data(stream & idx3==j).name)...
                    'Record(s):\n' ...
                    sprintf(' %d ',find(stream & idx3==j))]);
                badcmp=true;
            end
            if(numel(unique(cmpinc(stream & idx3==j)))>1)
                warning('seizmo:findtriplets:cmpChangesAzimuth',...
                    ['Inclination for component not consistent!\n' ...
                    'Filename(s):\n' ...
                    sprintf(' %s\n',data(stream & idx3==j).name)...
                    'Record(s):\n' ...
                    sprintf(' %d ',find(stream & idx3==j))]);
                badcmp=true;
            end
        end
        if(badcmp)
            offset=offset+1;
            continue;
        end
        
        % check E/N cmp are horizontal
        if(any((strcmpi(cmpstr(stream),'N') ...
                | strcmpi(cmpstr(stream),'E')) & cmpinc(stream)~=90))
            warning('seizmo:findtriplets:northEastNotHorizontal',...
                ['North/East components not oriented horizontally!\n' ...
                'Filename(s):\n' sprintf(' %s\n',data(stream).name) ...
                'Record(s):\n' sprintf(' %d ',find(stream))]);
            offset=offset+1;
            continue;
        end
        
        % check Z cmp is vertical
        if(any(strcmpi(cmpstr(stream),'Z') ...
                & cmpinc(stream)~=0 & cmpinc(stream)~=180))
            warning('seizmo:findtriplets:zNotVertical',...
                ['Z component not oriented vertically!' ...
                'Filename(s):\n' sprintf(' %s\n',data(stream).name) ...
                'Record(s):\n' sprintf(' %d ',find(stream))]);
            offset=offset+1;
            continue;
        end
        
        % check orthogonality between each component
        [rows,cols]=find(tril(true(nstreamcmp),-1));
        cmporient=[cmpinc(stream) cmpaz(stream)];
        cmporient=cmporient(idx4,:);
        if(~all(isorthogonal(cmporient(rows,:),cmporient(cols,:))))
            warning('seizmo:findtriplets:nonOrthogonal',...
                ['Components for this triplet are not orthogonal!\n' ...
                'Filename(s):\n' sprintf(' %s\n',data(stream).name) ...
                'Record(s):\n' sprintf(' %d ',find(stream))]);
            offset=offset+1;
            continue;
        end
        
        % require exactly 3 cmp in stream (separate <3 & >3)
        if(nstreamcmp<3)
            warning('seizmo:findtriplets:tooFewCmp',...
                ['Less than 3 components for this triplet!\n' ...
                'Filename(s):\n' sprintf(' %s\n',data(stream).name) ...
                'Record(s):\n' sprintf(' %d ',find(stream))]);
            offset=offset+1;
            continue;
        elseif(nstreamcmp>3)
            warning('seizmo:findtriplets:tooManyCmp',...
                ['More than 3 components for this triplet!\n' ...
                'Filename(s):\n' sprintf(' %s\n',data(stream).name) ...
                'Record(s):\n' sprintf(' %d ',find(stream))]);
            offset=offset+1;
            continue;
        end
        
        % update triplet number & allow export
        idx2(stream)=i-offset;
        istriplet(stream)=true;
    end

    % return indices of triplets
    idx1=idx1(istriplet);
    idx2=idx2(istriplet);
    idx3=idx3(istriplet);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
