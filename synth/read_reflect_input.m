function [o]=read_reflect_input(file)
%READ_REFLECT_INPUT    Reads a reflectivity input file
%
%    Usage:    s=read_reflect_input(file)
%
%    Description:
%     S=READ_REFLECT_INPUT(FILE) reads a reflectivity input file formatted
%     equivalently to that produced by MAKE_REFLECT_INPUT.  The output is a
%     struct S that contains all the reflectivity input information in
%     several fields.  No format checking is done: it will either go
%     smoothly or bust.
%
%    Notes:
%     - THIS READS A LOCAL FORMAT!  Do not expect this to read reflectivity
%       input files for your reflectivity programs!
%     - The reflect code that works with these functions can be found in
%       the same directory as this mfile in the gzipped tarball
%       reflect.tar.gz.  Please don't expect me to get this running for
%       you.
%
%    Examples:
%
%    See also: MAKE_REFLECT_INPUT, READ_REFLECT_OUTPUT, REFLECT2SEIZMO

%     Version History:
%        Aug. 10, 2010 - initial version
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 23:00 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% get ascii lines
if(nargin<1); file=[]; end
lines=getwords(readtxt(file,...
    {'*.dat;*.DAT' 'Reflect Input Files (*.dat,*.DAT)';
    '*.*' 'All Files (*.*)'}),sprintf('\n'));

% parse lines
% first 3 lines are identifiers: run, event, model names
tmp=strtrim(lines(1:3));
[o.run,o.event,o.model]=deal(tmp{:});
% nlayers, nstf, nsta, field, response, space, output, function, radiation
tmp=num2cell(str2double(cellstr(reshape(lines{4},4,[])')));
[o.nlayers,o.nsrcfun,o.nsta,o.field,o.response,o.space,o.synth,...
    o.function,o.radiation]=deal(tmp{:});
% npts, log2npts, delta, expfac
tmp=num2cell(str2double(cellstr(reshape(lines{5},10,4)')));
[o.npts,o.log2npts,o.delta,o.expfac]=deal(tmp{:});
% cmtlayer, fcmtdep
tmp=num2cell(str2double(cellstr(reshape(lines{6},10,2)')));
[o.cmtlayer,o.fcentroiddep]=deal(tmp{:});
% moment tensor
o.momten=reshape(str2double([cellstr(reshape(lines{7},10,3)');
                             cellstr(reshape(lines{8},10,3)');
                             cellstr(reshape(lines{9},10,3)')]),3,3)';
% freqlimits
o.freqlimits=str2double(cellstr(reshape(lines{10},10,2)'))';
% nslow, slowlimits
tmp=str2double(cellstr(reshape(lines{11},10,3)'));
o.nslow=tmp(1); o.slowlimits=[tmp(2) tmp(3)];
% denseslowflag, denseslowlimits
tmp=str2double(cellstr(reshape(lines{12},10,3)'));
o.dsflag=tmp(1); o.dslimits=[tmp(2) tmp(3)];
% slowtaperflag, slowtaperlimits
tmp=str2double(cellstr(reshape(lines{13},10,3)'));
o.taperflag=tmp(1); o.taperlimits=[tmp(2) tmp(3)];
% model - reflectivity, vp, vs, rho, thick, 1/qp, 1/qs
a=14+o.nlayers;
tmp=char(lines(14:a-1));
o.modval=...
    [str2double(cellstr(tmp(:,1))) str2double(cellstr(tmp(:,2:10))) ...
    str2double(cellstr(tmp(:,11:20))) str2double(cellstr(tmp(:,21:30))) ...
    str2double(cellstr(tmp(:,31:40))) str2double(cellstr(tmp(:,41:50))) ...
    str2double(cellstr(tmp(:,51:60)))];
% vazflag, filter, units
tmp=num2cell(str2double(cellstr(reshape(lines{a},4,3)')));
[o.vazflag,o.filter,o.units]=deal(tmp{:});
% stations
b=a+o.nsta;
tmp=char(lines(a+1:b));
o.staname=strtrim(cellstr(tmp(:,1:20)));
o.dist=str2double(cellstr(tmp(:,21:30)));
o.az=str2double(cellstr(tmp(:,31:40)));
o.stla=str2double(cellstr(tmp(:,41:50)));
o.stlo=str2double(cellstr(tmp(:,51:60)));
% srcfun (imaginary component dropped!)
tmp=str2double(cellstr(reshape(lines{b+1},10,[])'));
o.srcfun=tmp(1:2:end);
% vr, vrstart
tmp=num2cell(str2double(cellstr(reshape(lines{b+2},10,2)')));
[o.vr,o.vrstart]=deal(tmp{:});

end
