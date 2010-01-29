function [data]=prep_corrections(indir,outdir)

% get event dirs in indir
dirs=xdir(indir);
dirs={dirs.name}';
dirs(strcmp(dirs,'.') | strcmp(dirs,'..'))=[];

% loop over dirs
for i=1:numel(dirs)
	disp(dirs{i})
	data=readseizmo([indir filesep dirs{i}]);
	% add arrivals
	data=addarrivals(data,'phases','Pdiff,P','model','prem');
	% add corrections
	data=addcorrection(data,['corrections' filesep dirs{i} filesep 'Pdiff.mantle.hmsl-p06.upswing'],'resp0');
	data=addcorrection(data,['corrections' filesep dirs{i} filesep 'Pdiff.prem.crust'],'resp1');
	data=addcorrection(data,['corrections' filesep dirs{i} filesep 'Pdiff.prem.elev'],'resp2');
	data=addcorrection(data,['corrections' filesep dirs{i} filesep 'Pdiff.ellip.ak135'],'resp3');
	data=syncrates(data,5);
	writeseizmo(data,'path',[outdir filesep dirs{i}]);
end

end

