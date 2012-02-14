function [cmp]=version_compare(ver1,ver2)
%VERSION_COMPARE    Compares versions strings given in XX.XX.XX.... format
%
%    Usage:    cmp=version_compare(ver1,ver2)
%
%    Description:
%     CMP=VERSION_COMPARE(VER1,VER2) compares the versions in VER1 & VER2
%     returning CMP=1 if VER1>VER2, CMP=-1 if VER1<VER2, or CMP=0 if
%     VER1=VER2.  VER1 & VER2 must be strings.  Versions are split based on
%     the '.' delimiter (typically versions are formatted as major.minor
%     and so on).  Alphabetical characters have higher values than numbers
%     so that 'a'>'9'.
%
%    Notes:
%
%    Examples:
%     % Compare some simple versions:
%     cmp=version_compare('1','3')
%
%     % A case from wgrib2:
%     cmp=version_compare('v0.1.5f','v0.1.9.4')
%
%    See also: VERLESSTHAN, VER, PARSE_ALPHANUMERIC

%     Version History:
%        Feb. 13, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 13, 2012 at 13:00 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% require strings
if(~(ischar(ver1) || iscellstr(ver1)) ...
        || ~(ischar(ver2) || iscellstr(ver2)))
    error('seizmo:version_compare:badInput',...
        'VER must be a string!');
end

% now convert to cell strings
ver1=cellstr(ver1);
ver2=cellstr(ver2);

% expand scalars
if(isscalar(ver1)); ver1=ver1(ones(numel(ver2),1),1); end
if(isscalar(ver2)); ver2=ver2(ones(numel(ver1),1),1); end
ncmp=numel(ver1);

% loop over each row in ver1/ver2
cmp=nan(ncmp,1);
for i=1:ncmp
    % require single row
    if(size(ver1{i},1)>1 || size(ver2{i},1)>1)
        error('seizmo:version_compare:badInput',...
            'VER string must be a row vector!');
    end
    
    % get major.minor.revision...
    f1=getwords(ver1{i},'.');
    f2=getwords(ver2{i},'.');
    
    % loop over comparible fields until there is a difference
    nf1=numel(f1); nf2=numel(f2);
    for j=1:min(nf1,nf2)
        % try parsing as numbers
        d1=str2double(f1{j});
        d2=str2double(f2{j});
        
        if(~isnan(d1) && ~isnan(d2))
            % both are numbers
            cmp(i)=sign(d1-d2);
        else
            % parse as alphanumeric
            [an1,isnum1]=parse_alphanumeric(f1{j});
            [an2,isnum2]=parse_alphanumeric(f2{j});
            nbit1=numel(isnum1);
            nbit2=numel(isnum2);
            
            % loop over parsed bits
            for k=1:min(nbit1,nbit2)
                % n vs n
                if(isnum1(k) && isnum2(k))
                    cmp(i)=sign(an1{k}-an2{k});
                % a vs a
                elseif(~isnum1(k) && ~isnum2(k))
                    last=max(numel(an1{k}),numel(an2{k}));
                    for l=1:last
                        cmp(i)=sign(an1{k}(l)-an2{k}(l));
                        if(cmp(i)); break; end
                    end
                    % no diff so use number of alpha
                    if(~cmp(i))
                        cmp(i)=sign(numel(an1{k})-numel(an2{k}));
                    end
                % n vs a
                elseif(isnum1(k) && ~isnum2(k))
                    % alpha wins (0-9, a-z)
                    cmp(i)=-1;
                % a vs n
                elseif(~isnum1(k) && isnum2(k))
                    % alpha wins (0-9, a-z)
                    cmp(i)=1;
                end
                if(cmp(i)); break; end
            end
        end
        
        % continue until 1 or -1
        if(cmp(i)); break; end
    end
    
    % failed to find a difference in comparible fields
    % so just use the number of fields
    if(isnan(cmp(i)) || ~cmp(i)); cmp(i)=sign(nf1-nf2); end
end

end
