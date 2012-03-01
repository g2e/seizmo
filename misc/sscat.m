function [s]=sscat(varargin)
%SSCAT    Concatenates struct(s) into a scalar struct
%
%    Usage:    s=sscat(s1,s2,...)
%
%    Description:
%     S=SSCAT(S1,S2,...) concatenates structs S1, S2, etc and combines
%     their elements to create a scalar struct S.  Inputs are required to
%     be struct, have the same fields, and all fields should have the same
%     number of rows.
%
%    Notes:
%
%    Examples:
%     % Combine the results of a few CMT searches:
%     cmt1=findcmt('time',[2000 1]);
%     cmt2=findcmt('depth',700);
%     cmt3=findcmts('mwrange',[7 8]);
%     cmt=sscat(cmt1,cmt2,cmt3);
%
%    See also: SSIDX, READNDK, PARSE_ISC_ORIGIN, READSODEVENTCSV

%     Version History:
%        Feb. 29, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 29, 2012 at 17:25 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check inputs
for i=1:nargin
    if(~isstruct(varargin{i}))
        error('seizmo:sscat:badInput',...
            'SSCAT inputs must be struct!');
    end
    if(i==1)
        f=fieldnames(varargin{i});
    elseif(~isempty(setxor(fieldnames(varargin{i}),f)))
        error('seizmo:sscat:badInput',...
            'All structs must have the same fields!');
    end
end

% concatenate
try
    for i=1:nargin
        for j=1:numel(varargin{i})
            if(i==1 && j==1)
                for k=1:numel(f)
                    s.(f{k})=varargin{i}(j).(f{k});
                end
            else
                for k=1:numel(f)
                    s.(f{k})=[s.(f{k}); varargin{i}(j).(f{k})];
                end
            end
        end
    end
catch
    warning('seizmo:sscat:badInput',...
        'Failed to concatenate inputs!');
    error(lasterror);
end

end
