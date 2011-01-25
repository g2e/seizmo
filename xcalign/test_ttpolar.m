function [varargout]=test_ttpolar(p,w,n)
%TEST_TTPOLAR    Tests the success rate of the TTPOLAR algorithm
%
%    Usage:    test_ttpolar
%              test_ttpolar(p,w,n)
%              s=test_ttpolar(...)
%
%    Description:
%     TEST_TTPOLAR tests the TTPOLAR algorithm success rate and plots the
%     results.  It tests a range of % defects in a polarity grid (from 0%
%     to 50%) for a 50x50 polarity grid.  The test is repeated 100 times
%     for each % defective point.
%
%     TEST_TTPOLAR(P,W,N) allows changing the % defect steps (P), the range
%     of polarity grid sizes (W), or the number of tests per percent defect
%     point (N).  The defaults are P=0:.25:50, W=50, N=100.
%
%     S=TEST_TTPOLAR(...) returns S, a Nx2 matrix of [%defect %success].
%     The first column %defect is equal to the input P if given.
%
%    Notes:
%
%    Examples:
%     % Check out the relative success rate between several matrix sizes:
%     s10=test_ttpolar([],10);
%     s20=test_ttpolar([],20);
%     s30=test_ttpolar([],30);
%     s40=test_ttpolar([],40);
%     s50=test_ttpolar([],50);
%     figure;
%     plot(s10(:,1),s10(:,2),'rx:',...
%          s20(:,1),s20(:,2),'gp-.',...
%          s30(:,1),s30(:,2),'bv-',...
%          s40(:,1),s40(:,2),'ms--',...
%          s50(:,1),s50(:,2),'ko:');
%     xlabel('% Defects');
%     ylabel('% Success');
%     title('TTPOLAR Algorithm Performance');
%     legend({'10x10' '20x20' '30x30' '40x40' '50x50'});
%
%    See also: TTPOLAR

%     Version History:
%        Jan. 19, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 19, 2011 at 10:35 GMT

% todo

% check nargin
error(nargchk(0,3,nargin));

% defaults
if(nargin<1 || isempty(p)); p=0:.25:50; end
if(nargin<2 || isempty(w)); w=[50 50]; end
if(nargin<3 || isempty(n)); n=100; end

% check inputs
if(~isreal(p) || any(p<0 | p>100))
    error('seizmo:test_ttpolar:badInput',...
        'P must be percentages between 0 & 100!');
elseif(~isreal(w) || ~any(numel(w)==[1 2]) ...
        || any(w<=0) || any(w~=fix(w)))
    error('seizmo:test_ttpolar:badInput',...
        ['W must be a 1x2 vector giving the range ' ...
        'of allowed square matrix sizes!']);
elseif(~isreal(n) || ~isscalar(n) || n<=0 || n~=fix(n))
    error('seizmo:test_ttpolar:badInput',...
        'N must be the number of tests per percentage point!');
end

% expand scalar w
if(isscalar(w)); w=[w w]; end

% verbose messaging
npp=numel(p);
verbose=seizmoverbose;
if(verbose)
    disp('Testing TTPOLAR Algorithm');
    print_time_left(0,npp);
end

% loop over defect amounts
results=zeros(npp,1);
for a=1:npp
    for b=1:n
        % create randomly sized polarity grid
        sol=sign(rand(floor(0.5+w(1)+rand*diff(w)),1)-.5);
        pg=sol*sol.';
        
        % defect locations
        np=size(pg,1);             % number of rows
        np2=(np^2-np)/2;           % number of lower triangle elements
        lti=randperm(np2);         % lower triange indices
        [r,c]=lti2sub(np,lti);     % row/column subscripts
        z=round(np2*p(a)/100);     % number of defective positions
        idx1=sub2ind([np np],r(1:z),c(1:z)); % linear indices
        idx2=sub2ind([np np],c(1:z),r(1:z)); % linear indices
        
        % insert defects
        pg(idx1)=pg(idx1)*-1;
        pg(idx2)=pg(idx2)*-1;
        
        % solve
        pol=ttpolar(pg);
        
        % check
        results(a)=results(a)+(isequal(pol,sol) || isequal(pol,-1*sol));
    end
    
    % update counter
    if(verbose); print_time_left(a,npp); end
end
results=results/n*100;

% plot results
if(nargout)
    varargout{1}=[p(:) results];
else
    figure;
    plot(p,results,'kx:');
    xlabel('% Defects');
    ylabel('% Success');
    title('TTPOLAR Algorithm Performance');
end

end
