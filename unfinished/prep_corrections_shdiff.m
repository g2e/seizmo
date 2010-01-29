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
	data=addarrivals(data,'phases','Sdiff,S','model','prem');
	% add corrections
	data=addcorrection(data,['corrections' filesep dirs{i} filesep 'SHdiff.mantle.hmsl-s06.upswing'],'resp0');
	data=addcorrection(data,['corrections' filesep dirs{i} filesep 'Sdiff.prem.crust'],'resp1');
	data=addcorrection(data,['corrections' filesep dirs{i} filesep 'Sdiff.prem.elev'],'resp2');
	data=addcorrection(data,['corrections' filesep dirs{i} filesep 'SHdiff.ellip.ak135'],'resp3');
	data=syncrates(data,5);
	writeseizmo(data,'path',[outdir filesep dirs{i}]);
end

end

