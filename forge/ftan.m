function []=ftan(data)

% basically something like this:
adata=idft(omegaanalytic(omegagaussian(dft(rtr(data(100))),1./(5:.5:50),10)),'nonsymmetric');
imagesc(records2mat(cut(solofun(adata,@abs),-500,500)))

alpha=100;
bank=filter_bank([0.01 0.1],'variable',0.2,0.1)
w=bank(:,1);
fdata=dft(data);
z=nan(numel(x),numel(w));
for i=1:numel(w)
    fdata=omegagaussian(fdata,w,alpha);
    hdata=omegahilbert(fdata);
    f=ifft(fdata);
    h=ifft(hdata);
    
end

end
