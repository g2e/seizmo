function []=ftan(data)

% basically something like this:
delta=gh(data,'delta'); fnyq=1/(2*delta);
adata=idft(omegaanalytic(omegagaussian(dft(rtr(data)),linspace(fnyq/250,fnyq,249),100)),'nonsymmetric');
figure; imagesc(records2mat(solofun(adata,@abs)));

%{
function [x,fc]=ftan(x,fc,a,fs)
error(nargchk(1,4,nargin));
if(nargin<2 || isempty(fc)); fc=100; end
if(nargin<3 || isempty(a)); a=100; end
if(nargin<4 || isempty(fs)); fs=1; end
if(isscalar(fc)); fc=fs/2*(1/fc:1/fc:1); end
x=fft(x);

fdata=dft(data);
z=nan(numel(x),numel(w));
for i=1:numel(w)
    fdata=omegagaussian(fdata,w,alpha);
    hdata=omegahilbert(fdata);
    f=ifft(fdata);
    h=ifft(hdata);
    
end
%}
end
