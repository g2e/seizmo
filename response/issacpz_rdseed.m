function [lgc]=issacpz_rdseed(pz,type)
%ISSACPZ_RDSEED    Returns TRUE if input is a SAC PoleZero RDSEED struct
%
%    Usage:    lgc=issacpz_rdseed(pz)
%              lgc=issacpz_rdseed(pz,type)
%
%    Description:
%     LGC=ISSACPZ_RDSEED(PZ) checks for valid fieldnames and data types as
%     returned by READSACPZ_RDSEED.  Currently this assumes that the
%     comment block corresponds to that output by RDSEED v5.1.
%
%     LGC=ISSACPZ_RDSEED(PZ,TYPE) allows verifying a SAC PoleZero struct as
%     output by a different version of RDSEED or as the old database format
%     used by SEIZMO.  TYPE may currently be one of the following:
%      '5.1' - RDSEED 5.1 comment block (default)
%      'OLD' - SEIZMO SACPZDB struct format
%
%    Notes:
%     RDSEED v5.1 SACPZ data structure layout (scalar struct):
%      path             - path to PoleZero file
%      name             - PoleZero file name
%      knetwk           - network associated with PoleZero file
%      kstnm            - station associated with PoleZero file
%      kcmpnm           - component associated with PoleZero file
%      khole            - stream associated with PoleZero file
%      created          - time PoleZero file was created
%      b                - begin time for PoleZero validity
%      e                - end time for PoleZero validity
%      description      - station location description
%      stla             - station latitude in degrees
%      stlo             - station longitude in degrees
%      stel             - station elevation in meters
%      stdp             - station depth in meters
%      cmpinc           - component inclination in degrees from vertical
%      cmpaz            - component azimuth in degrees from North
%      sr               - sample rate in Hz
%      input            - input data units
%      output           - output data units
%      insttype         - type of instrument
%      instgain         - gain factor of instrument
%      instgainunits    - units associated with instrument gain factor
%      comment          - addition instrument comments
%      sensitivity      - total sensitivity (multiplication of gains)
%      sensitivityunits - units associated with sensitivity
%      a0               - PoleZero normalization factor (gives 1 at Fref)
%      z                - zeros in angular frequency
%      p                - poles in angular frequency
%      k                - constant (a0 * sensitivity)
%
%     Old SACPZ data structure layout (one pz per struct index):
%      path   - path to PoleZero file
%      name   - PoleZero file name
%      knetwk - network associated with PoleZero file
%      kstnm  - station associated with PoleZero file
%      kcmpnm - component associated with PoleZero file
%      khole  - stream associated with PoleZero file
%      b      - begin time for PoleZero validity
%      e      - end time for PoleZero validity
%      z      - zeros
%      p      - poles
%      k      - constant
%
%    Examples:
%     % Is the latest output from IRIS web services compatible?
%     url=['http://service.iris.edu/irisws/sacpz/1/' ...
%          'query?net=IU&loc=*&cha=*&sta=ANMO'];
%     pz=readsacpz_rdseed(urlread(url),true);
%     issacpz_rdseed(pz)
%
%    See also: READSACPZ_RDSEED, WRITESACPZ_RDSEED, READSACPZ, WRITESACPZ,
%              REMOVESACPZ, APPLYSACPZ, MAKESACPZDB, GENSACPZNAME,
%              PARSE_SACPZ_FILENAME, GETSACPZ, FIX_OLD_SACPZ

%     Version History:
%        Feb. 25, 2014 - initial version
%        Mar.  5, 2014 - 5.3 changed to 5.1
%        Mar.  6, 2014 - update See also section
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  6, 2014 at 15:05 GMT

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
valid.TYPE={'5.1' 'OLD'};
if(nargin<2 || isempty(type)); type='5.1'; end
if(~isstring(type) || ~any(strcmpi(type,valid.TYPE)))
    error('seizmo:issacpz_rdseed:badInput',...
        'TYPE is not a valid string!');
end

% requirements by type
switch upper(type)
    case '5.1'
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
