function [phaseout,repetitions]=mkstriprepetitions(phase);
% mkstriprepetitions.......determine number of repetitions in multiple phase
%
% call: [phaseout,repetitions]=mkstriprepetitions(phase);
%
%       phase: string containing a seismic phase name
%
% result: phaseout: as the input PHASE, but with all trailing digits removed
%
%         repetitions: the trailing digits that were initially contained in
%                      PHASE.
%
%
% For the handling of multiple phases Xn, where X is any seismic phase name
% and n is the "degree" of the multiple, it is necessary to determine the
% number n and the phase on which the multiple is based.
% If, for example the seismic phase X is SKS, then SKS2 denotes the second
% multiple SKSSKS. MKSTRIPREPETITIONS returns the string "SKS" as base
% phase and the number 2 as number of repetitions.
%
% If PHASE does not end with a digit, 0 is returned as number of
% repetitions.
%
% Limitations:
%
% The phase PK5P2 denotes PKKKKKPPKKKKKP, according to the IASPEI
% nomenclature. However, this routine does not resolve the "5" and returns
% "PK5P" as base phase and 2 as the number of repetitions.
%
% Martin Knapmeyer, 20.02.2006


%%% init result
phaseout=phase;
repetitions=[];

%%% transform input into ascii codes
phasesav=phase;
phase=abs(phase);

if (phase(end)>=48)&(phase(end)<=57)
    %%% PHASE ends with a number - just a quick identification to save time
    repstring='';
    done=0;
    len=length(phase);
    indy=len;
    while done==0
        if (phase(indy)>=48)&(phase(indy)<=57)
            repstring=strvcat(char(phase(indy)),repstring);
            indy=indy-1;
            if indy==0
                done=1;
            end; % if indy==0
        else
            done=1;
        end; % if (phase(indy)>=48)&(phase(indy)<=57)
    end; % while done
    repetitions=str2num(repstring');
    phaseout=phaseout(1:(indy));
else
    repetitons=0;
end; % if (phase(end)>=48)&(phase(end)<=57)
