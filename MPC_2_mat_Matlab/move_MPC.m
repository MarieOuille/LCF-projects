maxw=max(max(abs(uw.^2)));
maxt=max(max(abs(u.^2)));

figure(5);

dim = [0.2 0.2 0.2 0.2];

for kk=1:50:nz
    
    an=annotation('textbox',dim,'Str',['z = ' num2str(z(kk)*1e2,'%.1f') ' cm'],'FitBoxToText','on');

    subplot(211)
    %plot(lambdav*1e9,10*log10(abs(uw(kk,:).^2)./maxw));
    plot(lambdav(lambdav>0)*1e9,abs(uw(kk,lambdav>0).^2)./maxw);
    xlabel('Wavelength (nm)');
    ylabel('Spectrum (a. u.)');
    axis([3000 7000 0 1.1]);
    %axis([800 1500 -100 0]);
    
    subplot(212)
    plot(t*1e15,abs(u(kk,:).^2)./maxt);
    xlabel('Time (fs)');
    ylabel('Intensity (a. u.)');
    axis([-1000 1000 0 1.1]);
    
    drawnow
    delete(an)

end
