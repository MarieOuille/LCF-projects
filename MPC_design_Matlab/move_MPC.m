maxw=max(max(abs(uw.^2)));
maxt=max(max(abs(u.^2)));

figure(5);

for kk=1:10:nz
    
    subplot(211)
    %plot(lambdav*1e9,10*log10(abs(uw(kk,:).^2)./maxw));
    plot(lambdav(lambdav>0)*1e9,abs(uw(kk,lambdav>0).^2)./maxw,'k');
    xlabel('Wavelength (nm)');
    ylabel('Spectrum (a. u.)');
    axis([1000 2500 0 1.1]);
    %axis([800 1500 -100 0]);
    
    subplot(212)
    plot(t*1e15,abs(u(kk,:).^2)./maxt,'k');
    xlabel('Time (fs)');
    ylabel('Intensity (a. u.)');
    axis([-500 500 0 1.1]);
    
    drawnow

end
