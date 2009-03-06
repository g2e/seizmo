function [f,amp]=getpeak(freq,totwin,delta,window,tprlen)
% makes a sensitivity kernel for two-plane wave method

% generate sin function
t=totwin(1):delta:totwin(2);
x=sin(2*pi*freq*t);

% window indices
idx=(t>=window(1) & t<=window(2));
nidx=sum(idx);

% cosine taper
tprr=0.5*cos(2*pi*linspace(0,0.5,tprlen/delta))+0.5;
tpr=tprr(end:-1:1);
if(numel(tpr)*2>nidx); error('window too small'); end
coswin=[tpr ones(1,nidx-numel(tpr)*2) tprr];

% take window
x(~idx)=0;
x(idx)=x(idx).*coswin;

% take fft
npts=numel(x);
nfft=2^nextpow2(npts);
f=1/(delta*2)*linspace(0,1,nfft/2);
s=fft(x,nfft)/npts;

% get amplitude
amp=2*abs(s(1:nfft/2));

% pull central peak
[dumb,idx]=max(amp);
damp=diff(amp);
idx1=find(damp(1:idx-1)<0,1,'last')+1;
idx2=idx+find(damp(idx+1:end)>0,1,'first');
if(isempty(idx1)); idx1=1; end
if(isempty(idx2)); idx2=numel(amp); end
f=f(idx1:idx2);
amp=amp(idx1:idx2);

end