function [x,t]=make_source_timefunction(delta,hwidth,type)
%MAKE_SOURCE_TIMEFUNCTION    Returns source time functions for convolution
%
%    Usage: [x,t]=make_source_timefunction(delta,hwidth)
%           [x,t]=make_source_timefunction(delta,hwidth,type)
%
%    Description: [X,T]=MAKE_SOURCE_TIMEFUNCTION(DELTA,HWIDTH) creates
%     gaussian source functions with a halfwidth of HWIDTH sampled at a
%     spacing of DELTA.  DELTA must be a real array with values >0.  HWIDTH
%     must be a real array.  DELTA and HWIDTH should be equal sized or be
%     scalar.  Scalar expansion is enabled.  The gaussian values are
%     returned in X and the times are in T.  X and T are cell arrays with
%     the source function corresponding to each element in DELTA/HWIDTH.
%
%     [X,T]=MAKE_SOURCE_TIMEFUNCTION(DELTA,HWIDTH,TYPE) sets the type of
%     the source function.  Valid values of TYPE are 'GAUSSIAN' &
%     'TRIANGLE'.  Type may be a cell array of strings.
%
%    Notes:
%     - gaussian-type functions extend from -1.5*HWIDTH to 1.5*HWIDTH
%     - triangle-type functions extend from -HWIDTH to HWIDTH
%
%    Examples:
%     Make source functions for all records in a SEIZMO dataset:
%      delta=getheader(data,'delta');
%      [x,t]=make_source_timefunction(delta,10);
%
%     Now plot all the functions together:
%      tmp=[t x]';
%      plot(tmp{:});
%
%    See also: CONVOLVE_SOURCE_TIMEFUNCTION, TRIANGLETF, GAUSSIANTF,
%              DECONVOLVE_SOURCE_TIMEFUNCTION

%     Version History:
%        Oct. 17, 2009 - initial version
%        Oct. 29, 2009 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 29, 2009 at 01:50 GMT

% todo:

% check nargin
msg=nargchk(2,3,nargin);
if(~isempty(msg)); error(msg); end

% check delta
if(isempty(delta) || ~isreal(delta) || any(delta<=0))
    error('seizmo:make_source_timefunction:badInput',...
        'DELTA must be real and >0!');
else
    delta=delta(:);
    nd=numel(delta);
end

% check hwidth
if(isempty(hwidth) || ~isreal(hwidth))
    error('seizmo:make_source_timefunction:badInput',...
        'HWIDTH must be real!');
else
    hwidth=hwidth(:);
    nh=numel(hwidth);
end

% check type
valid={'TRIANGLE' 'GAUSSIAN'};
if(nargin<3 || isempty(type)); type={'GAUSSIAN'}; end
if(ischar(type)); type=cellstr(type); end
if(~iscellstr(type) || any(~ismember(upper(type),valid)))
    error('seizmo:make_source_timefunction:badInput',...
        ['TYPE must be a char/cellstr array containing:\n' ...
        '''GAUSSIAN'' or ''TRIANGLE''!']);
else
    type=type(:);
    nt=numel(type);
end

% check/expand
nstf=max([nd nt nh]);
if(~all([nd nt nh]==1 | [nd nt nh]==nstf))
    error('seizmo:make_source_timefunction:badInput',...
        ['DELTA, TYPE, & HWIDTH must be scalar\n' ...
        'or have an equal number of elements!']);
end
if(nd==1); delta=delta(ones(nstf,1),1); end
if(nt==1); type=type(ones(nstf,1),1); end
if(nh==1); hwidth=hwidth(ones(nstf,1),1); end

% gaussian parameters
a=exp(0.5)*sqrt(2);
b=exp(0.5)/sqrt(pi);

% loop over each
x=cell(nstf,1);
t=cell(nstf,1);
for i=1:nstf
    switch upper(type{i})
        case 'GAUSSIAN'
            npts=ceil(1.5*hwidth(i)/delta(i));
            t{i}=delta(i).*(-npts:npts);
            x{i}=gaussiantf(t{i},0,hwidth(i),a,b/hwidth(i));
        case 'TRIANGLE'
            npts=ceil(hwidth(i)/delta(i));
            t{i}=delta(i).*(-npts:npts);
            x{i}=triangletf(t{i},0,hwidth(i),1/hwidth(i));
    end
end

end
