function [sod]=cmt2sod(cmt,type,mag)
%CMT2SOD    Returns a SOD event csv struct from a GlobalCMT struct
%
%    Usage:    sod=cmt2sod(cmt)
%              sod=cmt2sod(cmt,type)
%              sod=cmt2sod(cmt,type,mag)
%
%    Description:
%     SOD=CMT2SOD(CMT) extracts the hypocenter info for a GlobalCMT
%     centroid moment tensor and puts it into a SOD event struct.  This is
%     useful for handing off to WRITESODEVENTCSV to facilitate programs
%     that need this info in that format.
%
%     SOD=CMT2SOD(CMT,TYPE) specifies if the output event info is based on
%     the hypocenter or the centroid.  TYPE may be 'HYPO' or 'CENTROID'.
%     The default is 'HYPO'.
%
%     SOD=CMT2SOD(CMT,TYPE,MAG) specifies the magnitude type to output.
%     The default is 'MW' but 'MS' or 'MB' are also valid.
%
%    Notes:
%
%    Examples:
%     % Create a SOD event csv file for all available GlobalCMT CMTs:
%     writesodeventcsv('globalcmt.csv',cmt2sod(findcmts));
%
%    See also: SOD2CMT, FINDCMTS, FINDCMT,
%              READSODEVENTCSV, WRITESODEVENTCSV

%     Version History:
%        Feb. 29, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 29, 2012 at 17:25 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% defaults
if(nargin<2 || isempty(type)); type='hypo'; end
if(nargin<3 || isempty(mag)); mag='mw'; end

% check inputs
if(~isstruct(cmt) && ~isscalar(cmt))
    error('seizmo:cmt2sod:badInput',...
        'CMT must be a scalar struct formatted like FINDCMTS output!');
elseif(~ischar(type) || ndims(type)>2 || size(type,1)>1)
    error('seizmo:cmt2sod:badInput',...
        'TYPE must be ''HYPO'' or ''CENTROID''!');
elseif(~ischar(mag) || ndims(mag)>2 || size(mag,1)>1)
    error('seizmo:cmt2sod:badInput',...
        'MAG must be ''MB'' ''MS'' or ''MW''!');
elseif(~any(strcmpi(mag,{'MB' 'MS' 'MW'})))
    error('seizmo:cmt2sod:badInput',...
        'MAG must be ''MB'' ''MS'' or ''MW''!');
end

% process depending on type/mag
switch lower(type)
    case {'h' 'hypo' 'hypocenter'}
        % hypo
        sod.time=[cmt.year cmt.month cmt.day cmt.hour ...
            cmt.minute cmt.seconds];
        sod.latitude=cmt.latitude;
        sod.longitude=cmt.longitude;
        sod.depth=cmt.depth;
        sod.depthUnits={'kiloMETER'};
        sod.depthUnits=sod.depthUnits(ones(size(sod.depth)));
        switch lower(mag)
            case 'mb'
                sod.magnitude=cmt.mb;
                sod.magnitudeType={'MB'};
            case 'ms'
                sod.magnitude=cmt.ms;
                sod.magnitudeType={'MS'};
            case 'mw'
                sod.magnitude=momentmag(cmt);
                sod.magnitudeType={'MW'};
        end
        sod.magnitudeType=sod.magnitudeType(ones(size(sod.magnitude)));
        sod.catalog=cmt.catalog;
    case {'c' 'cmt' 'centroid'}
        sod.time=fixtimes([cmt.year cmt.month cmt.day cmt.hour ...
            cmt.minute cmt.seconds+cmt.centroidtime]);
        sod.latitude=cmt.centroidlat;
        sod.longitude=cmt.centroidlon;
        sod.depth=cmt.centroiddep;
        sod.depthUnits={'kiloMETER'};
        sod.depthUnits=sod.depthUnits(ones(size(sod.depth)));
        switch lower(mag)
            case 'mb'
                sod.magnitude=cmt.mb;
                sod.magnitudeType={'MB'};
            case 'ms'
                sod.magnitude=cmt.ms;
                sod.magnitudeType={'MS'};
            case 'mw'
                sod.magnitude=momentmag(cmt);
                sod.magnitudeType={'MW'};
        end
        sod.magnitudeType=sod.magnitudeType(ones(size(sod.magnitude)));
        sod.catalog={'CMT'};
        sod.catalog=sod.catalog(ones(size(sod.depth)));
    otherwise
        error('seizmo:cmt2sod:badInput',...
            'TYPE Unrecognized: %s',type);
end

end
