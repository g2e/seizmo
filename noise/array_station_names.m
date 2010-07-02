function [names]=array_station_names(array)
%ARRAY_STATION_NAMES    Returns station names for an array
%
%    Usage:    arrays=array_station_names()
%              stations=array_station_names(array)
%
%    Description: ARRAYS=ARRAY_STATION_NAMES() returns the available array
%     names for this function (you can add more if you want - email me and
%     I will perma-add them).
%
%     STATIONS=ARRAY_STATION_NAMES(ARRAY) returns station names for the
%     array given by ARRAY.  ARRAY does not have to be and is very likely
%     not a network code.  This is useful for extracting specific station
%     sets (subarrays & virtual arrays) from data.
%
%    Notes:
%     - Subarrays are definitely supported.  For example, I have split the
%       Ethiopia array XI (2000-2002) that was actually 2 arrays into their
%       Ethiopia and Kenya portions.
%
%    Examples:
%     What arrays are available?
%      array_station_names
%
%    See also: DAYDIRS_STACKCORR

%     Version History:
%        June 27, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 27, 2010 at 02:25 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% available array names
arrays={'cameroon' 'ethiopia' 'kenya' 'south africa' 'tanzania'};

% if no inputs return arrays
if(~nargin)
    % display or return as cell array
    if(~nargout)
        fprintf('Available Arrays:\n');
        fprintf('%s\n',arrays{:});
    else
        names=arrays;
    end
    return;
end

% check array name
if(~ischar(array) || ndims(array)~=2 || size(array,1)~=1)
    error('seizmo:array_station_names:badInput',...
        'ARRAY must be a string!');
end

% arrays
cameroon={'CM01' 'CM02' 'CM03' 'CM04' 'CM05' 'CM06' 'CM07' 'CM08' ...
          'CM09' 'CM10' 'CM11' 'CM12' 'CM13' 'CM14' 'CM15' 'CM16' ...
          'CM17' 'CM18' 'CM19' 'CM20' 'CM21' 'CM22' 'CM23' 'CM24' ...
          'CM25' 'CM26' 'CM27' 'CM28' 'CM29' 'CM30' 'CM31' 'CM32'};
ethiopia={'AAUS' 'ARBA' 'BAHI' 'BDAR' 'BELA' 'BIRH' 'BUTA' ...
          'CHEF' 'DELE' 'DIYA' 'DMRK' 'FICH' 'GOBA' 'GUDE' ...
          'HERO' 'HIRN' 'HOSA' 'JIMA' 'KARA' 'NAZA' 'NEKE' ...
          'SELA' 'TEND' 'TERC' 'WANE' 'WASH' 'WELI' 'WELK'};
kenya={'ANGA' 'BARI' 'BOKO' 'KAKA' 'KIBO' ...
       'KITU' 'KR42' 'NARO' 'NDEI' 'TALE'};
south_africa={'SA00' 'SA01' 'SA02' 'SA03' 'SA04' 'SA05' 'SA07' ...
              'SA08' 'SA09' 'SA10' 'SA11' 'SA12' 'SA13' 'SA139' ...
              'SA14' 'SA15' 'SA155' 'SA16' 'SA169' 'SA17' 'SA18' ...
              'SA19' 'SA20' 'SA22' 'SA23' 'SA24' 'SA25' 'SA26' ...
              'SA27' 'SA28' 'SA29' 'SA30' 'SA31' 'SA32' 'SA33' ...
              'SA34' 'SA35' 'SA36' 'SA37' 'SA38' 'SA39' 'SA40' ...
              'SA42' 'SA43' 'SA44' 'SA45' 'SA46' 'SA47' 'SA48' ...
              'SA49' 'SA50' 'SA51' 'SA52' 'SA53' 'SA54' 'SA55' ...
              'SA56' 'SA57' 'SA58' 'SA59' 'SA60' 'SA61' 'SA62' ...
              'SA63' 'SA64' 'SA65' 'SA66' 'SA67' 'SA68' 'SA69' ...
              'SA70' 'SA71' 'SA72' 'SA73' 'SA74' 'SA75' 'SA76' ...
              'SA77' 'SA78' 'SA79' 'SA80' 'SA81' 'SA82'};
tanzania={'AMBA' 'BASO' 'GOMA' 'HALE' 'INZA' 'KIBA' 'KIBE' ...
          'KOMO' 'KOND' 'LONG' 'MBWE' 'MITU' 'MTAN' 'MTOR' ...
          'PAND' 'PUGE' 'RUNG' 'SING' 'TARA' 'TUND' 'URAM'};

switch lower(array)
    case 'cameroon'
        names=cameroon;
    case 'ethiopia'
        names=ethiopia;
    case 'kenya'
        names=kenya;
    case 'south africa'
        names=south_africa;
    case 'tanzania'
        names=tanzania;
    otherwise
        error('seizmo:arraynames:badname',...
            'ARRAY %s UNKNOWN!',upper(array));
end

end
