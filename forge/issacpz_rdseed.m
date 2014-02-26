function [lgc]=issacpz_rdseed(pz,type)
%ISSACPZ_RDSEED    Returns TRUE if input is a SAC PoleZero RDSEED struct
%
%    Usage:    lgc=issacpz_rdseed(pz)
%              lgc=issacpz_rdseed(pz,type)
%
%    Description:
%     LGC=ISSACPZ_RDSEED(PZ) checks for valid fieldnames and data types as
%     returned by READSACPZ_RDSEED.  Currently this assumes that the
%     comment block corresponds to that output by RDSEED v5.3.
%
%     LGC=ISSACPZ_RDSEED(PZ,TYPE) allows verifying a SAC PoleZero struct as
%     output by a different version of RDSEED or as the old database format
%     used by SEIZMO.  TYPE may currently be one of the following:
%      '5.3' - RDSEED 5.3 comment block (default)
%      'OLD' - SEIZMO SACPZDB struct format
%
%    Notes:
%
%    Examples:
%     % Is the latest output from IRIS web services compatible?
%     url=['http://service.iris.edu/irisws/sacpz/1/' ...
%          'query?net=IU&loc=*&cha=*&sta=ANMO'];
%     pz=readsacpz_rdseed(urlread(url),true);
%     issacpz_rdseed(pz)
%
%    See also: READSACPZ_RDSEED, WRITESACPZ_RDSEED, READSACPZ, WRITESACPZ,
%              REMOVESACPZ, APPLYSACPZ, MAKESACPZDB, DB2SACPZ,
%              GENSACPZNAME, PARSE_SACPZ_FILENAME, GETSACPZ, FIX_OLD_SACPZ

%     Version History:
%        Feb. 25, 2014 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 25, 2014 at 15:05 GMT

% todo:

% check number of inputs
error(nargchk(1,2,nargin));

% require pz is a struct
lgc=true;
if(~isstruct(pz))
    lgc=false;
    return;
end

% default/check type
valid.TYPE={'5.3' 'OLD'};
if(nargin<2 || isempty(type)); type='5.3'; end
if(~isstring(type) || ~any(strcmpi(type,valid.TYPE)))
    error('seizmo:issacpz_rdseed:badInput',...
        'TYPE is not a valid string!');
end

% requirements by type
switch upper(type)
    case '5.3'
        fields={'path' 'name' 'knetwk' 'kstnm' 'khole' 'kcmpnm' ...
            'created' 'b' 'e' 'description' 'stla' 'stlo' 'stel' 'stdp' ...
            'cmpinc' 'cmpaz' 'sr' 'input' 'output' 'insttype' ...
            'instgain' 'instgainunits' 'comment' 'sensitivity' ...
            'sensitivityunits' 'a0' 'z' 'p' 'k'};
        floats=[0 0 0 0 0 0 1 1 1 0 1 1 1 1 1 1 1 0 0 0 1 0 0 1 0 1 1 1 1];
        celled=[1 1 1 1 1 1 0 0 0 1 0 0 0 0 0 0 0 1 1 1 0 1 1 0 1 0 1 1 0];
    case 'OLD'
        fields={'path' 'name' 'knetwk' 'kstnm' 'khole' 'kcmpnm' ...
            'b' 'e' 'z' 'p' 'k'};
        floats=[0 0 0 0 0 0 1 1 1 1 1];
        celled=[0 0 0 0 0 0 0 0 0 0 0];
end

% test requirements
if(~all(ismember(fields,fieldnames(pz))))
    lgc=false;
    return;
end
for i=1:numel(floats)
    if(floats(i))
        if(celled(i))
            if(~iscell(pz(1).(fields{i})))
                lgc=false;
                return;
            elseif(~isfloat(pz(1).(fields{i}){1}))
                lgc=false;
                return;
            end
        else
            if(~isfloat(pz(1).(fields{i})))
                lgc=false;
                return;
            end
        end
    else
        if(celled(i))
            if(~iscellstr(pz(1).(fields{i})))
                lgc=false;
                return;
            end
        else
            if(~isstring(pz(1).(fields{i})))
                lgc=false;
                return;
            end
        end
    end
end

end
