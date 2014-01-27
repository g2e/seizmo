function []=write_1dmodel_nd(file,model,overwrite,noheader)
%WRITE_1DMODEL_ND    Writes a 1D model struct in .nd format
%
%    Usage:    write_1dmodel_nd(file,model)
%              write_1dmodel_nd(file,model,overwrite)
%              write_1dmodel_nd(file,model,overwrite,noheader)
%
%    Description:
%     WRITE_1DMODEL_ND(FILE,MODEL) writes an .nd formatted file (see
%     references for the origins of the format) using the information in
%     the struct MODEL.  The .nd format stands for "named discontinuity"
%     and as the name implies, the major discontinuities (moho, cmb, icb)
%     are named (mantle, outer-core, inner-core).  Two comment lines are
%     inserted at the top of the model file and contain the model name and
%     some of the flags (ocean, crust, isotropic, refperiod, flattened).
%     FILE is a string giving the path and name of the output .nd file.
%     MODEL must follow the struct layout as output from PREM, IASP91,
%     AK135, etc.  The model must be isotropic!
%
%     WRITE_1DMODEL_ND(FILE,MODEL,OVERWRITE) quietly overwrites the pre-
%     existing ND file without confirmation when OVERWRITE is set to TRUE.
%     By default OVERWRITE is FALSE.  OVERWRITE is ignored in the graphical
%     file creation menu.
%
%     WRITE_1DMODEL_ND(FILE,MODEL,OVERWRITE,NOHEADER) skips writing the
%     header portion of the .nd file when NOHEADER is set to TRUE.  This is
%     helpful for programs like TAUP that cannot handle comments.  By
%     default NOHEADER is FALSE.
%
%    Notes:
%     - References for the .nd file format:
%        Davis, J. P. and I. H. Henson (1993a). Development of an X-Windows
%         tool to compute Gaussian bean synthetic seismograms. Technical
%         Report TGAL-93-03, Phillip Laboratory, Hancom AFB, MA.
%        Davis, J. P. and I. H. Henson (1993b). User’s Guide to Xgbm: An
%         X-Windows System to compute Gaussian bean synthetic seismograms
%         (1.1 ed.). Alexandria, VA: Teledyne Geotech Alexandria
%         Laboratories.
%
%    Examples:
%     % Write out PREM:
%     write_1dmodel_nd('prem_test.nd',prem);
%
%    See also: PREM, IASP91, AK135, PREM_PERFECT, PREM2_PERFECT,
%              PERTURB_1DMODEL, PLOT1DMODEL, CMB_1DMODEL_LIBRARY
%              FLATTEN_1DMODEL, CHK1DMODEL, MODELS

%     Version History:
%        Sep. 18, 2010 - initial version
%        Sep. 19, 2010 - support for inf Q output as 0
%        Feb. 21, 2012 - noheader flag
%        Jan. 26, 2014 - abs path exist fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2014 at 10:35 GMT

% todo

% check nargin
error(nargchk(2,4,nargin));

% directory separator
fs=filesep;

% default overwrite to false
if(nargin<3 || isempty(overwrite)); overwrite=false; end
if(nargin<4 || isempty(noheader)); noheader=false; end
if(~isscalar(overwrite) || ~islogical(overwrite))
    error('seizmo:write_1dmodel_nd:badInput',...
        'OVERWRITE flag must be a scalar logical!');
elseif(~isscalar(noheader) || ~islogical(noheader))
    error('seizmo:write_1dmodel_nd:badInput',...
        'NOHEADER flag must be a scalar logical!');
end

% check model
error(chk1dmodel(model));

% has to be isotropic
if(~model.isotropic)
    error('seizmo:write_1dmodel_nd:badInput',...
        'MODEL must be isotropic for .nd format!');
end

% graphical selection
if(isempty(file))
    [file,path]=uiputfile(...
        {'*.ND;*.nd' 'Named-Discontinuity Files (*.ND,*.nd)';
        '*.*' 'All Files (*.*)'},...
        'Save 1D model as');
    if(isequal(0,file))
        error('seizmo:write_1dmodel_nd:noFileSelected',...
            'No output file selected!');
    end
    file=[path fs file];
else
    % check file
    if(~isstring(file))
        error('seizmo:write_1dmodel_nd:fileNotString',...
            'FILE must be a string!');
    end
    if(~isabspath(file)); file=[pwd fs file]; end
    if(exist(file,'file'))
        if(exist(file,'dir'))
            error('seizmo:write_1dmodel_nd:dirConflict',...
                '1D Model .ND File: %s\nIs A Directory!',file);
        end
        if(~overwrite)
            fprintf('1D Model .ND File: %s\nFile Exists!\n',file);
            reply=input('Overwrite? Y/N [N]: ','s');
            if(isempty(reply) || ~strncmpi(reply,'y',1))
                disp('Not overwriting!');
                return;
            end
            disp('Overwriting!');
        end
    end
end

% get discontinuities
dctopidx=find(diff(model.depth)==0);
dcdepth=model.depth(dctopidx);

% find major discontinuity depths
% moho == mantle
% cmb  == outer-core
% icb  == inner-core
% moho is first discontinuity <100km where lowerside Vp >= 7.4km/s
% cmb is discontinuity 2870-2910km where lowerside Vs <= 1.0km/s
% icb is discontinuity 5100-5200km where lowerside Vs >= 1.0km/s
moho=dctopidx(dcdepth<100 ...
    & model.vp(dctopidx)<=7.4 & model.vp(dctopidx+1)>=7.4);
cmb=dctopidx(dcdepth>2870 & dcdepth<2910 ...
    & model.vs(dctopidx)>0 & model.vs(dctopidx+1)==0);
icb=dctopidx(dcdepth>5130 & dcdepth<5170 ...
    & model.vs(dctopidx)==0 & model.vs(dctopidx+1)>0);

% fail if a major discontinuity is not found
if(~any(moho) || ~any(cmb) || ~any(icb))
    error('seizmo:write_1dmodel_nd:badInput',...
        'Could not locate all major discontinuities!');
end

% open file for writing as ascii
fid=fopen(file,'wt');

% check if file is openable
if(fid<0)
    error('seizmo:write_1dmodel_nd:cannotOpenFile',...
        '1D Model .ND File: %s\nNot Openable!',file);
end

% write 2 comment lines with metadata
if(~noheader)
    fprintf(fid,'# %s\n',model.name);
    fprintf(fid,'# %d %d %d %d %d\n',...
        model.ocean,model.crust,model.isotropic,...
        model.refperiod,model.flattened);
end

% write model
fields={'depth' 'vp' 'vs' 'rho' 'qk' 'qu' 'qp' 'qs'};
props=ismember(fields,fieldnames(model));
if(isequal(props,[1 1 1 0 0 0 0 0]))
    % just depth vp vs
    fprintf(fid,'%9.3f %10.5f %10.5f\n',...
        [model.depth(1:moho) model.vp(1:moho) model.vs(1:moho)]');
    fprintf(fid,'mantle\n');
    fprintf(fid,'%9.3f %10.5f %10.5f\n',[model.depth(moho+1:cmb) ...
        model.vp(moho+1:cmb) model.vs(moho+1:cmb)]');
    fprintf(fid,'outer-core\n');
    fprintf(fid,'%9.3f %10.5f %10.5f\n',[model.depth(cmb+1:icb) ...
        model.vp(cmb+1:icb) model.vs(cmb+1:icb)]');
    fprintf(fid,'inner-core\n');
    fprintf(fid,'%9.3f %10.5f %10.5f\n',[model.depth(icb+1:end) ...
        model.vp(icb+1:end) model.vs(icb+1:end)]');
elseif(isequal(props,[1 1 1 1 0 0 0 0]))
    % just depth vp vs rho
    fprintf(fid,'%9.3f %10.5f %10.5f %10.5f\n',...
        [model.depth(1:moho) model.vp(1:moho) ...
        model.vs(1:moho) model.rho(1:moho)]');
    fprintf(fid,'mantle\n');
    fprintf(fid,'%9.3f %10.5f %10.5f %10.5f\n',...
        [model.depth(moho+1:cmb) model.vp(moho+1:cmb) ...
        model.vs(moho+1:cmb) model.rho(moho+1:cmb)]');
    fprintf(fid,'outer-core\n');
    fprintf(fid,'%9.3f %10.5f %10.5f %10.5f\n',...
        [model.depth(cmb+1:icb) model.vp(cmb+1:icb) ...
        model.vs(cmb+1:icb) model.rho(cmb+1:icb)]');
    fprintf(fid,'inner-core\n');
    fprintf(fid,'%9.3f %10.5f %10.5f %10.5f\n',...
        [model.depth(icb+1:end) model.vp(icb+1:end) ...
        model.vs(icb+1:end) model.rho(icb+1:end)]');
elseif(isequal(props,[1 1 1 1 1 1 0 0]))
    % depth vp vs rho qk qu
    % convert qk qu vp vs to qp
    qp=qkqu2qp(model.qk,model.qu,model.vp,model.vs);
    % change inf qp qs to 0
    qp(isinf(qp))=0;
    model.qu(isinf(model.qu))=0;
    fprintf(fid,'%9.3f %10.5f %10.5f %10.5f %10.1f %10.1f\n',...
        [model.depth(1:moho) model.vp(1:moho) model.vs(1:moho) ...
        model.rho(1:moho) qp(1:moho) model.qu(1:moho)]');
    fprintf(fid,'mantle\n');
    fprintf(fid,'%9.3f %10.5f %10.5f %10.5f %10.1f %10.1f\n',...
        [model.depth(moho+1:cmb) model.vp(moho+1:cmb) ...
        model.vs(moho+1:cmb) model.rho(moho+1:cmb) ...
        qp(moho+1:cmb) model.qu(moho+1:cmb)]');
    fprintf(fid,'outer-core\n');
    fprintf(fid,'%9.3f %10.5f %10.5f %10.5f %10.1f %10.1f\n',...
        [model.depth(cmb+1:icb) model.vp(cmb+1:icb) ...
        model.vs(cmb+1:icb) model.rho(cmb+1:icb) ...
        qp(cmb+1:icb) model.qu(cmb+1:icb)]');
    fprintf(fid,'inner-core\n');
    fprintf(fid,'%9.3f %10.5f %10.5f %10.5f %10.1f %10.1f\n',...
        [model.depth(icb+1:end) model.vp(icb+1:end) ...
        model.vs(icb+1:end) model.rho(icb+1:end) ...
        qp(icb+1:end) model.qu(icb+1:end)]');
elseif(isequal(props,[1 1 1 1 1 1 1 1]) ...
        || isequal(props,[1 1 1 1 0 0 1 1]))
    % depth vp vs rho qp qs
    % change inf qp qs to 0
    model.qp(isinf(model.qp))=0;
    model.qs(isinf(model.qs))=0;
    fprintf(fid,'%9.3f %10.5f %10.5f %10.5f %10.1f %10.1f\n',...
        [model.depth(1:moho) model.vp(1:moho) model.vs(1:moho) ...
        model.rho(1:moho) model.qp(1:moho) model.qs(1:moho)]');
    fprintf(fid,'mantle\n');
    fprintf(fid,'%9.3f %10.5f %10.5f %10.5f %10.1f %10.1f\n',...
        [model.depth(moho+1:cmb) model.vp(moho+1:cmb) ...
        model.vs(moho+1:cmb) model.rho(moho+1:cmb) ...
        model.qp(moho+1:cmb) model.qs(moho+1:cmb)]');
    fprintf(fid,'outer-core\n');
    fprintf(fid,'%9.3f %10.5f %10.5f %10.5f %10.1f %10.1f\n',...
        [model.depth(cmb+1:icb) model.vp(cmb+1:icb) ...
        model.vs(cmb+1:icb) model.rho(cmb+1:icb) ...
        model.qp(cmb+1:icb) model.qs(cmb+1:icb)]');
    fprintf(fid,'inner-core\n');
    fprintf(fid,'%9.3f %10.5f %10.5f %10.5f %10.1f %10.1f\n',...
        [model.depth(icb+1:end) model.vp(icb+1:end) ...
        model.vs(icb+1:end) model.rho(icb+1:end) ...
        model.qp(icb+1:end) model.qs(icb+1:end)]');
else
    error('seizmo:write_1dmodel_nd:unsupported',...
        'Unsupported set of model properties for .nd file creation!');
end

% close file
fclose(fid);

end
