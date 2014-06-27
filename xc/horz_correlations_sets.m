function [in,set,cmp,rev]=horz_correlations_sets(xc)
%HORZ_CORRELATIONS_SETS   Returns indices for horiz. correlation sets
%
%    Usage:    [in,set,cmp]=horz_correlations_sets(xc)
%              [in,set,cmp,rev]=horz_correlations_sets(xc)
%
%    Description:
%     [IN,SET,CMP]=HORZ_CORRELATIONS_SETS(XC) groups horiz. correlograms
%     into rotatible sets for ROTATE_CORRELATIONS.  This is similar to the
%     functions HORZPAIRS & FINDTRIPLETS which find the rotatible data for
%     the functions ROTATE/ROTATE3.  XC is a SEIZMO struct of correlograms
%     created by CORRELATE.  Sets are identified by having common master &
%     slave stations, horizontal orientation, orthogonality, and lag times
%     (lag times must be equal based on the B, E, NPTS & DELTA fields while
%     for frequency domain data this is also based on SB, SDELTA, & NSPTS -
%     therefore time domain and freqeuncy domain data are never grouped).
%     Sets may be of size 3 (only for auto-correlation of a single stations
%     components) or 4.  IN gives the indices of the correlograms in XC
%     that are in a set.  SET gives the corresponding set indices (use
%     max(SET) to get the number of sets). CMP gives the corresponding
%     component indices (1, 2, 3, or 4).  Component indices are as follows:
%          (1) NN or RR
%          (2) NE or RT
%          (3) EN or TR
%          (4) EE or TT
%     Auto-correlation sets may potentially be missing component 2 or 3 due
%     to the redundant information (2=reverse(3) or reverse(2)=3).
%
%     [IN,SET,CMP,REV]=HORZ_CORRELATIONS_SETS(XC) also indicates what
%     correlograms that are in the horizontal sets need to be reversed.
%     This is useful when master & slave stations are switched for a set.
%     REV is a logical vector equal in size to IN where each value
%     indicates whether the corresponding record given in IN needs to be
%     reversed.  For example, if REV(3)=TRUE then XC(IN(3)) should be
%     reversed using REVERSE_CORRELATIONS.
%
%    Notes:
%     - Actually doesn't require E/N or R/T, just that the underlying
%       components are orthogonal.
%     - Set separation based on orientation will detect false positives for
%       sets that share the same station pair, lags and have equal or
%       orthogonal orientation.  In this case the sets will be ignored and
%       not included in the output.  So don't throw a bunch of NE & RT sets
%       for the same station pairs together and expect it all to be
%       separated cleanly!  This also means you should not expect this to
%       parse a full matrix of horizontal correlations as that will have 8
%       correlations per pair (4 of which are redundant) - please use
%       NO_REDUNDANT_CORRELATIONS beforehand!
%
%    Examples:
%     % Remove vertical correlations and non-set horizontal correlations:
%     xc=xc(horz_correlations_sets(xc));
%
%     % Perform reversal to make the sets match:
%     [in,set,cmp,rev]=horz_correlations_sets(xc);
%     xc(in(rev))=reverse_correlations(xc(in(rev)));
%
%     % Get sets and loop over them:
%     for i=1:max(set)
%         % record indices for this pair
%         ridx=in(set==i);
%
%         % separate indices based on component
%         cidx1=in(set==i & cmp==1);
%         cidx2=in(set==i & cmp==2);
%         cidx3=in(set==i & cmp==3);
%         cidx4=in(set==i & cmp==4);
%
%         ... % do something useful here % ...
%     end
%
%    See also: CORRELATE, SPLIT_AUTO_CORRELATIONS, ROTATE_CORRELATIONS,
%              REVERSE_CORRELATIONS, ISXC, NO_REDUNDANT_CORRELATIONS,
%              NAME_CORRELATIONS, IS_FULL_MATRIX_OF_CORRELATIONS,
%              HORZPAIRS, FINDTRIPLETS

%     Version History:
%        Mar.  9, 2013 - initial concept drafted
%        Sep.  5, 2013 - first working version
%        Sep. 20, 2013 - properly optimized checking, deg fudge factor
%        Jan. 15, 2014 - fixed several typos of this functions name
%        June  4, 2014 - doc update, position field checks
%        June 12, 2014 - handle fd i/o
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check structure & headers
xc=checkheader(xc,...
    'MULCMP_DEP','ERROR',...
    'XYZ_IFTYPE','ERROR',...
    'FALSE_LEVEN','ERROR',...
    'UNSET_ST_LATLON','ERROR',...
    'UNSET_EV_LATLON','ERROR');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% safely get necessary header info
try
    [kuser,b,e,npts,delta,snm,mnm,...
        scmp,mcmp,st,ev,delaz,sb,sdelta,nspts]=getheader(xc,...
        'kuser','b','e','npts','delta','kname','kt','cmp','user',...
        'st','ev','delaz','sb','sdelta','nspts');
    mcmp=mcmp(:,3:4); % drop the record indices from correlate
    se=sb+sdelta.*(nspts-1);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% number of records
nrecs=numel(xc);

% error if non-correlations
if(~all(strcmp(kuser(:,1),'MASTER') & strcmp(kuser(:,2),'SLAVE')))
    error('seizmo:horz_correlations_sets:badInput',...
        'Some records appear not to be correlations!');
end

% only horizontal correlations
horz=mcmp(:,1)==90 & scmp(:,1)==90;
if(sum(horz)<3); [in,set,cmp,rev]=deal([]); return; end

% default output
% - in will be converted from logical to linear and the other
%   arrays will be subsetted to the corresponding elements
in=false(nrecs,1);
set=nan(nrecs,1);
cmp=nan(nrecs,1);
rev=false(nrecs,1);

% stream names
mnm=lower(strcat(mnm(:,1),'.',mnm(:,2),'.',...
    mnm(:,3),'.',strnlen(char(mnm(:,4)),2)));
snm=lower(strcat(snm(:,1),'.',snm(:,2),'.',...
    snm(:,3),'.',strnlen(char(snm(:,4)),2)));

% autocorrelations
auto=strcmp(mnm,snm);
if(any(auto) && any(b(auto)~=-e(auto)))
    error('seizmo:horz_correlations_sets:autoxcWithAsymetricLags',...
        'Autocorrelation sets must not have asymmetric lags!');
end

% division into supersets based on naming & lags
% - 2 superset names to handle master<=>slave ambiguity
ssetname1=strcat(mnm,'_',snm,'_',...
    num2str(sb),'_',num2str(nspts),'_',num2str(sdelta),'_',...
    num2str(b),'_',num2str(npts),'_',num2str(delta));
ssetname2=strcat(snm,'_',mnm,'_',...
    num2str(-se),'_',num2str(nspts),'_',num2str(sdelta),'_',...
    num2str(-e),'_',num2str(npts),'_',num2str(delta));
[ssetname,r2ssidx,ss2ridx]=unique([ssetname1;ssetname2],'first');

% number of supersets
nssets=numel(ssetname);

% double over set to record indexing so we
% know the 2 supersets each record is in
ss2ridx=reshape(ss2ridx,nrecs,2);

% loop over each record and find sets
% - removing redundant supersets due to master<=>slave ambiguity
% - spliting supersets into sets based on cmpaz
ckdss=false(nssets,1); % logical indicating if the superset is checked
nsets=0; % number of sets
for i=1:nrecs
    % current superset
    css=ss2ridx(i,1);
    
    % already checked this superset?
    if(ckdss(css)); continue; end
    
    % mark this superset & its sister as checked
    % - this removes redundant supersets
    ckdss(ss2ridx(i,:))=true;
    
    % who is in this superset?
    inss=find((ss2ridx(:,1)==css | ss2ridx(:,2)==css) & horz);
    rev(ss2ridx(:,1)~=css & ss2ridx(:,2)==css)=true;
    sspop=numel(inss);
    
    % skip if not enough members to form valid set
    if(sspop<3); continue; end
    
    % form sets based on orientation
    inaset=false(sspop,1);
    for j=1:sspop
        % skip if already in a set
        if(inaset(j)); continue; end
        
        % require equal position fields for sanity
        ssst=st(inss,:);
        ssev=ev(inss,:);
        ssda=delaz(inss,:);
        if(any(rev(inss)))
            ssst(rev(inss),:)=ssev(rev(inss),:);
            ssev(rev(inss),:)=ssst(rev(inss),:);
            ssda(rev(inss),[3 2])=ssda(rev(inss),[2 3]);
        end
        if(size(unique(ssst,'rows'),1)~=1 ...
                || size(unique(ssev,'rows'),1)~=1)
            continue;
        end
        f=1e-8; % degrees fudge factor
        if(any(abs(ssda(:,1)-ssda(1,1))>f ...
                | abs(azdiff(ssda(:,2),ssda(1,2)))>f ...
                | abs(azdiff(ssda(:,3),ssda(1,3)))>f ...
                | abs(ssda(:,4)-ssda(1,4))>f))
            continue;
        end
        
        % grab azimuths (for ease)
        maz=mcmp(inss,2);
        saz=scmp(inss,2);
        if(any(rev(inss)))
            maz(rev(inss))=scmp(inss(rev(inss)),2);
            saz(rev(inss))=mcmp(inss(rev(inss)),2);
        end
        
        % who is in this set?
        % - based on orientation (parallel or orthogonal)
        % - note that this fails to detect if we have both +/-90
        %   - this would occur for multiple sets for this station pair that
        %     have component orientations that are orthogonal to the other
        %     sets orientation
        f=1e-8; % degrees fudge factor
        ins=(abs(azdiff(maz(j),maz))<f ...
            | abs(abs(azdiff(maz(j),maz))-90)<f) ...
            & (abs(azdiff(saz(j),saz))<f ...
            | abs(abs(azdiff(saz(j),saz))-90)<f);
        
        % mark members
        inaset(ins)=true;
        
        % autoxc set?
        as=auto(inss(find(ins,1,'first')));
        
        % too few/many correlograms in this set?
        setpop=sum(ins);
        if(as)
            % auto-correlation set requires 3 or 4
            if(setpop<3 || setpop>4); continue; end
        else
            % cross correlation set requires exactly 4
            if(setpop~=4); continue; end
        end
        
        % how many azimuths for each station?
        % - need exactly 2 for both
        % - this catches the orthogonal sets problem mentioned above
        %   - these are skipped so you probably want to avoid multiple
        %     sets of correlograms for the same station pair as this
        %     is a point of failure
        umaz=unique(maz(ins));
        usaz=unique(saz(ins));
        if(numel(umaz)~=2 || numel(usaz)~=2); continue; end
        
        % set up component indexing based on orientation convention
        if(abs(azdiff(umaz(1),umaz(2))-90)<f)
            % umaz(1) is "N/R" while umaz(2) is "E/T"
            mazi=[0 2];
        else % -90
            % umaz(1) is "E/T" while umaz(2) is "N/R"
            mazi=[2 0];
        end
        if(abs(azdiff(usaz(1),usaz(2))-90)<f)
            % usaz(1) is "N/R" while usaz(2) is "E/T"
            sazi=[1 2];
        else % -90
            % usaz(1) is "E/T" while usaz(2) is "N/R"
            sazi=[2 1];
        end
        
        % who is what component index-wise?
        cmpi=mazi((maz(ins)==umaz(2))+1)+sazi((saz(ins)==usaz(2))+1);
        
        % look out for redundant members
        if(numel(cmpi)~=numel(unique(cmpi))); continue; end
        
        % auto can only be missing cmp 2 or 3
        if(as && setpop==3)
            if(~any(cmpi==1) || ~any(cmpi==4)); continue; end
        end
        
        % okay set confirmed, update output
        nsets=nsets+1;
        ins=inss(ins); % logical to linear
        in(ins)=true;
        set(ins)=nsets;
        cmp(ins)=cmpi;
    end
end

% convert output
in=find(in); % logical to linear
set=set(in);
cmp=cmp(in);
rev=rev(in);

end
