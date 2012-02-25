function [tauobj]=taupcreate(filein,fileout)
%TAUPCREATE    Convert models to TauP toolkit format
%
%    Usage:    tauobj=taupcreate('filein','fileout')
%              tauobj=taupcreate('filein')
%              tauobj=taupcreate(1dmodel,'fileout')
%              tauobj=taupcreate(1dmodel)
%
%    Description:
%     TAUOBJ=TAUPCREATE('FILEIN','FILEOUT') reads the velocity model in
%     file FILEIN and converts it into TauP's .taup format, writing the
%     result as FILEOUT.  A java TauModel object is returned which can be
%     used for other TauP-based functions.
%
%     TAUOBJ=TAUPCREATE('FILEIN') skips writing the .taup file.
%
%     TAUOBJ=TAUPCREATE(1DMODEL,'FILEOUT') creates the TauP file/object
%     from a 1D model struct (eg. like that returned by PREM).
%
%     TAUOBJ=TAUPCREATE(1DMODEL) skips writing the .taup file.
%
%    Notes:
%
%    Examples:
%     % Create a SYLO model:
%     taupcreate(cmb_1dmodel_library('sylo'),'sylo.taup');
%
%    See also: TAUPTIME, TAUPPATH, TAUPCURVE, TAUPPIERCE, TAUP

%     Version History:
%        Feb. 24, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 24, 2012 at 10:35 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% try adding *.jar if no TauP class exists
if(~exist('edu.sc.seis.TauP.VelocityModel','class'))
    fs=filesep;
    mypath=fileparts(mfilename('fullpath'));
    javaaddpath([mypath fs 'lib' fs 'MatTauP-1.2beta4.jar']);
    javaaddpath([mypath fs 'lib' fs 'TauP-1.2beta4.jar']);
    javaaddpath([mypath fs 'lib' fs 'seisFile-1.0.8.jar']);
end

% get velocity model object
% - allows either a filename, vmobj, or 1dmodel struct
if(ischar(filein))
    % file on disk
    [filepath,filename,fileext]=fileparts(filein);
    if(~any(ismember(fileext,{'.nd' '.tvel'})))
        error('seizmo:taupcreate:badInput',...
            'Only .tvel or .nd files allowed!');
    end
    filein=javaMethod('readVelocityFile',...
        'edu.sc.seis.TauP.VelocityModel',...
        fullfile(filepath,[filename fileext]),fileext(2:end));
    filein.setModelName(filename); % should this be the 1st line?
    if(~filein.validate)
        error('seizmo:taupcreate:badInput',...
            'FILEIN velocity model not valid!');
    end
elseif(strcmp(class(filein),'edu.sc.seis.TauP.VelocityModel'))
    % velocity model object
    if(~filein.validate)
        error('seizmo:taupcreate:badInput',...
            'FILEIN velocity model object not valid!');
    end
elseif(isempty(chk1dmodel(filein)))
    % 1D model struct (write out & read in)
    newfile=tempname;
    modelname=filein.name;
    write_1dmodel_nd(newfile,filein,true,true);
    filein=javaMethod('readVelocityFile',...
        'edu.sc.seis.TauP.VelocityModel',newfile,'nd');
    filein.setModelName(modelname);
    if(~filein.validate)
        error('seizmo:taupcreate:badInput',...
            'FILEIN velocity model struct not valid!');
    end
    delete(newfile);
else
    error('seizmo:taupcreate:badInput',...
        'FILEIN must be a filename or a valid 1D model struct!');
end

% sample model
tcobj=javaObject('edu.sc.seis.TauP.TauP_Create');
tcobj.setVelocityModel(filein);
tauobj=tcobj.createTauModel;

% save model
if(nargin==2)
    if(ischar(fileout))
        tauobj.writeModel(fileout);
    else
        error('seizmo:taupcreate:badInput',...
            'FILEOUT must be a filename!');
    end
end

end
