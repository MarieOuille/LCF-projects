% Propagation in a capillary for radially symmetric modes
% Ionisation, instantaneous and delayed (Raman) Kerr, self-steepening,
% Vectorial
% Full waveguide and material dispersion
% Pressure gradient, mode resolved in matrix form
% Checked
% Adaptive step based on B integral

clearvars;
tic;

% Initial condition and capillary geometry
lambda0=1030e-9;
tgauss=28e-15/(sqrt(2*log(2))); 
energie=120e-6;                  
ell=1e4; % ellipticity 1e4 x 1e-4 y 1 circ

input_type='Gauss'; %Gauss or Exp_Speck or Exp_Ek
File='spectrum_measured_dscan_yoann.txt';

a=75e-6;
ncap=1.45;
pressurein=0.4; 
pressureout=0.4; 
gas='Ar';   
Z=1;

% Temporal and spatial ranges and steps

dt=0.2e-15;
T=1000e-15;

dz_min=1e-6;
dz_max=50e-3;
Bint_Kerr_max=0.1;
Bint_Ion_max=0.1;
nz_max=10000;

dr=1e-6; % high enough resolution to guarantee numerical mode orthogonality
nmode=5; % 1<nmode<20
SS_Ion_correction=1;

%Compression parameters
phi2negmax=-200e-30;
dphi2=5e-30;


%_________________________________________

% Parameters
c0=299792458;
charge=1.602e-19;
me=9.11e-31;
hbar=6.62e-34;
eps0=8.85e-12;
mu0=4*pi*1e-7;
w0=2*pi*c0/lambda0;
rgauss=0.64*a;
energie=energie*(1-exp(-2*a^2/rgauss^2)); % Gaussian coupling to capillary modes


% Time/frequency variables
%dt=1.01*lambda0/2/c0;
nt=2^(ceil(log2(T/dt)));
t=-(nt/2-1)*dt:dt:(nt/2)*dt;
w=2*pi*(1/dt)*(-nt/2:nt/2-1)./nt;
dw=w(2)-w(1);

lambdav=c0*2*pi./(w+w0);
lmicron=lambdav*1e6;
W=w+w0;

Nfield=2^4;
dtfield=dt/Nfield;
ntfield=nt*Nfield;
tfield=-(ntfield/2-1)*dtfield:dtfield:(ntfield/2)*dtfield;

% Transverse spatial variable
dr=a/round(a/dr);
r=0:dr:a;
nr=length(r);

% Capillary modes

[R,T]=ndgrid(r,t);
[Rfield,Tfield]=ndgrid(r,tfield);

u01=2.40482555769577;
u02=5.52007811028631;
u03=8.65372791291101;
u04=11.7915344390143;
u05=14.9309177084878;
u06=18.0710639679109;
u07=21.2116366298793;
u08=24.3524715307493;
u09=27.4934791320403;
u10=30.6346064684320;
u11=33.7758202135736;
u12=36.9170983536640;
u13=40.0584257646282;
u14=43.1997917131767;
u15=46.3411883716618;
u16=49.4826098973978;
u17=52.6240518411150;
u18=55.7655107550200;
u19=58.9069839260809;
u20=62.0484691902272;


ukm=[u01 u02 u03 u04 u05 u06 u07 u08 u09 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20];
ukm=ukm(1:nmode);
A=zeros(nmode,nr);
At=zeros(nmode,nr,nt);
alpha=zeros(nmode,nt);
Wm=ones(nmode,1)*w;

for km=1:nmode
    A(km,:)=besselj(0,ukm(km)*r/a);
    A(km,:)=A(km,:)./sqrt(sum(abs(A(km,:).^2).*2*pi.*r.*dr));
    alpha(km,:)=(ukm(km)/2/pi)^2*lambdav.^2/a^3*(0.5*(ncap^2+1)/sqrt(ncap^2-1));
    %alpha=0;
    At(km,:,:)=A(km,:).'*ones(1,nt);
end

alpha=real(alpha);

% Input field

switch input_type

    case 'Gauss'
    Ex=exp(-(R/rgauss).^2).*exp(-(T/tgauss).^2);
    energieint=sum(sum(abs(Ex.^2)*dt*2*pi.*R*dr));
    Ex=sqrt(energie/energieint).*Ex;

    case 'Exp_Ek'
    aa=load(File);
    tdata=aa(:,1)*1e-15;
    edata=sqrt(aa(:,2)).*exp(1i*aa(:,3));
    einterp=interp1(tdata,edata,t);
    einterp(isnan(einterp))=0;
    Ex=ones(nr,1)*einterp.*exp(-(R/rgauss).^2);
    energieint=sum(sum(abs(Ex.^2)*dt*2*pi.*R*dr));
    Ex=sqrt(energie/energieint).*Ex;
    
    case 'Exp_Speck'
    aaa=load(File);
    ldata=aaa(:,1)*1e-9;
    aaa(:,2) = aaa(:,2)./max(aaa(:,2));
    ewdata=sqrt(aaa(:,2)).*exp(-1i*aaa(:,3)); %might need to change phase sign
    ewdata=circshift(ewdata,0);
    ewinterp=interp1(ldata,ewdata,lambdav);
    ewinterp(isnan(ewinterp))=0;
    einterp=fftshift(ifft(ifftshift(ewinterp)));
    Ex=ones(nr,1)*einterp.*exp(-(R/rgauss).^2); 
    energieint=sum(sum(abs(Ex.^2)*dt*2*pi.*R*dr));
    Ex=sqrt(energie/energieint).*Ex;
end

Pcrete=max(squeeze(sum(abs(Ex.^2)*2*pi.*R*dr,1)));
Icrete=max(max(abs(Ex.^2)));

% Elliptical polarisation with axes aligned with x,y
Ex=ell/sqrt(1+ell^2)*Ex;
Ey=1i/ell*Ex;


ex=zeros(nmode,nt);
ey=zeros(nmode,nt);

for km=1:nmode
ex(km,:)=squeeze(sum(Ex.*conj(squeeze(At(km,:,:)))*2*pi.*R*dr,1));
ey(km,:)=squeeze(sum(Ey.*conj(squeeze(At(km,:,:)))*2*pi.*R*dr,1));
end

% Gas parameters

switch gas
    
    case 'Ar'
    nlambda=(2.50141e-3./(91.012-lmicron.^-2)+5.00283e-4./(87.892-lmicron.^-2)+5.22343e-2./(214.02-lmicron.^-2)); %Ar
    nlambda(isnan(nlambda))=0;
    nlambda=real(nlambda);
    n2_1bar=0.97e-23; 
    nstar=0.9289; 
    Ip=15.76*1.6e-19; 
    kappaio=9.18e16; 
    muio=1.28e12;
    C_TL=2.44;
    l_TL=1;
    Z_TL=1;
    alpha_TL=9.0;
    Ip_TL=Ip/4.35974417e-18; % in atomic units
    m_TL=0;
    kappa_TL=sqrt(2*Ip_TL);
    tau1=1e-15;
    tau2=1e-15;        
    fr=0;
    ion_lookup=load('linear_rate_Ar_m=0_lambda=1030nm.dat');
    %ion_lookup=load('linear_rate_Ar_m=0_lambda=515nm.dat');
    
    case 'Ne'
    nlambda=(0.00128145./(184.661-lmicron.^-2)+0.0220486./(376.840-lmicron.^-2));
    nlambda(isnan(nlambda))=0;
    nlambda=real(nlambda);
    n2_1bar=7.5e-25; 
    nstar=0.7942;
    Ip=21.56*1.6e-19;
    kappaio=1.14e17; 
    muio=2.05e12;
    C_TL=2.10;
    l_TL=1;
    Z_TL=1;
    alpha_TL=9.0;
    Ip_TL=Ip/4.35974417e-18; % in atomic units
    m_TL=0;
    kappa_TL=sqrt(2*Ip_TL);
    tau1=1e-15;
    tau2=1e-15;        
    fr=0;
    ion_lookup=load('linear_rate_Ne_m=0_lambda=1030nm.dat');
    
    case 'Xe'
    nlambda=(0.00322869./(46.301-lmicron.^-2)+0.00355393./(59.578-lmicron.^-2)+0.0606764./(112.74-lmicron.^-2)); %Xe
    nlambda(isnan(nlambda))=0;
    nlambda=real(nlambda);
    n2_1bar=5.2e-23; 
    nstar=1.0589;
    Ip=12.13*1.6e-19;
    kappaio=7.58e16;
    muio=8.65e11;
    C_TL=2.57;
    l_TL=1;
    Z_TL=1;
    alpha_TL=9.0;
    Ip_TL=Ip/4.35974417e-18; % in atomic units
    m_TL=0;
    kappa_TL=sqrt(2*Ip_TL);
    tau1=1e-15;
    tau2=1e-15;        
    fr=0;
    ion_lookup=load('linear_rate_Xe_m=0.dat');
    
    case 'He'
    nlambda=(0.014755297./(426.29740-lmicron.^-2));
    nlambda(isnan(nlambda))=0;
    nlambda=real(nlambda);
    n2_1bar=4.2e-25;
    nstar=0.7437;
    Ip=24.59*1.6e-19;
    kappaio=1.24e17;
    muio=2.50e12;
    C_TL=3.13;
    l_TL=0;
    Z_TL=1;
    alpha_TL=6.0;
    Ip_TL=Ip/4.35974417e-18; % in atomic units
    m_TL=0;
    kappa_TL=sqrt(2*Ip_TL);
    tau1=1e-15;
    tau2=1e-15;        
    fr=0;    
    ion_lookup=load('linear_rate_He_3_lambda=1030nm.dat');

    %%%%%%%%%%%%%%%%%%%%%%%%%
    case 'N2'
    nlambda=(6.8552e-5+3.243157e-2./(144-lmicron.^-2)); %N2
    nlambda(isnan(nlambda))=0;
    nlambda=real(nlambda);
    n2_1bar=4e-23; %N2 global
    tau1=66e-15; %62.5e-15 N2
    tau2=133e-15;  %120e-15 N2        
    fr=0.68;       %fr_1=0.75 for N2
    
    nstar=0.7437; %WARNING copied from He, must be found for N2
    Ip=24.59*1.6e-19;
    kappaio=1.24e17;
    muio=2.50e12;
    C_TL=3.13;
    l_TL=0;
    Z_TL=1;
    alpha_TL=6.0;
    Ip_TL=Ip/4.35974417e-18; % in atomic units
    m_TL=0;
    kappa_TL=sqrt(2*Ip_TL);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
end

Pcrit_1bar=lambdav(nt/2)^2/8/n2_1bar;
hr=(tau1^2+tau2^2)/(tau1*tau2^2).*exp(-t/tau2).*sin(t/tau1);
hr(t<0)=0;
hr=hr/sum(hr);
Resp=fr*hr;
Resp(t==0)=Resp(t==0)+(1-fr);
Respw=fft(Resp);
Respw=ones(nr,1)*Respw;

% Propagation

uxn=ex;
uyn=ey;
ux1=zeros(nz_max,nt);
uy1=zeros(nz_max,nt);

Emodex=zeros(nmode,nz_max);
Emodey=zeros(nmode,nz_max);
Euxtot=zeros(1,nz_max);
Euytot=zeros(1,nz_max);
Bintv_Kerr=zeros(nz_max,1);
Bintv_Ion=zeros(nz_max,1);
Ion_fraction=zeros(nz_max,1);

beta=zeros(nmode,nt);

Spectral_window=exp(-(w/(dw*nt*0.45)).^40);
Spectral_window=ones(nmode,1)*Spectral_window;

Temporal_window=exp(-(t/(dt*nt*0.45)).^40);
Temporal_window=ones(nmode,1)*Temporal_window;

z=zeros(1,nz_max);
z(1)=0;
index_z=1;
dz=dz_min;
FLAG=0;

while FLAG==0

waitbar(z(index_z)/Z);

pressure=sqrt(pressurein^2+z(index_z)/Z*(pressureout^2-pressurein^2));
ne0=2.688e25*pressure;
n=1+pressure*nlambda;
n2=pressure*n2_1bar;

for km=1:nmode
    beta(km,:)=W./c0.*sqrt(n.^2-ukm(km)^2*c0^2/a^2./W.^2);
end

beta=real(beta);
beta1=gradient(beta(1,:),W);
vg0=1/beta1(nt/2);



% Propagation over dz

   uxw=fftshift(fft(uxn,[],2),2);
   uyw=fftshift(fft(uyn,[],2),2);

   uxw=uxw.*Spectral_window;
   uyw=uyw.*Spectral_window;
   
   % Dispersion
   
   uxw=uxw.*exp(1i*dz/2.*(-beta+Wm/vg0)).*exp(-alpha*dz/2);
   uxw(:,W<0)=0;
   ux=ifft(ifftshift(uxw,2),[],2);
   
   uyw=uyw.*exp(1i*dz/2.*(-beta+Wm/vg0)).*exp(-alpha*dz/2);
   uyw(:,W<0)=0;
   uy=ifft(ifftshift(uyw,2),[],2);
   
   ux=ux.*Temporal_window;
   uy=uy.*Temporal_window;
   
   % Total field
  
   uxtot=zeros(nr,nt);
   uytot=zeros(nr,nt);
   
   for km=1:nmode
       uxtot=uxtot+ones(nr,1)*ux(km,:).*squeeze(At(km,:,:));
       uytot=uytot+ones(nr,1)*uy(km,:).*squeeze(At(km,:,:));
   end
   
   
   % Kerr + Raman nonlinearity
   
   zobiwanx=circshift(fftshift(ifft(fft(abs(uxtot.^2)+2/3*abs(uytot.^2),[],2).*Respw,[],2),2),1,2); %Kerr + Raman % add a coherent cross term for validity in the general polarization case ?
   zobiwany=circshift(fftshift(ifft(fft(abs(uytot.^2)+2/3*abs(uxtot.^2),[],2).*Respw,[],2),2),1,2); %Kerr + Raman
   
   uxtot=uxtot.*exp(-dz*1i*w0/c0*n2*zobiwanx);
   uytot=uytot.*exp(-dz*1i*w0/c0*n2*zobiwany);
   
   [Ftx,~]=gradient(+1i*w0/c0*n2*(zobiwanx.*uxtot+1/3*uytot.^2.*conj(uxtot)),dt,1);
   uxtot=uxtot+1i/w0*dz*Ftx-1i/3.*w0/c0*n2*uytot.^2.*conj(uxtot)*dz;
   
   [Fty,~]=gradient(+1i*w0/c0*n2*(zobiwany.*uytot+1/3*uxtot.^2.*conj(uytot)),dt,1);
   uytot=uytot+1i/w0*dz*Fty-1i/3.*w0/c0*n2*uxtot.^2.*conj(uytot)*dz; 
   
   % Ionization with cycle-averaged intensity lookup tables
   
%    efieldx=(interp1(t,uxtot.',tfield,'linear','extrap')).';
%    efieldy=(interp1(t,uytot.',tfield,'linear','extrap')).';
%    
%    efieldx=sqrt(2/eps0/c0)*real(efieldx.*exp(1i*w0*Tfield));
%    efieldy=sqrt(2/eps0/c0)*real(efieldy.*exp(1i*w0*Tfield));
%    
%    absefield=sqrt(efieldx.^2+efieldy.^2); %valeur absolue du champ en V/m
%    F_TL= absefield/5.14220652e11; % field strength (atomic units);
   
   % Wio2=kappaio*(muio./absefield).^(2*nstar-1).*exp(-muio/3./absefield); %these malvache pour comparaison
   
%    Wio2=C_TL^2 /2^abs(m_TL) /factorial(abs(m_TL))...% Tong Lin 2005 d'après HHGmax
%      * (2*l_TL+1)*factorial(l_TL+abs(m_TL)) /2 /factorial(l_TL-abs(m_TL))...
%      * 1 / kappa_TL^( 2*Z_TL/kappa_TL - 1 )...
%      * (2*kappa_TL^3 ./ F_TL) .^ ( 2*Z_TL/kappa_TL - abs(m_TL) - 1 )...
%     .* exp( -2/3 * kappa_TL^3./F_TL )...
%     .* exp( -alpha_TL * (Z_TL^2/Ip_TL) * (F_TL/kappa_TL^3) );
%    Wio2=Wio2*(1/2.418884326505e-17); %back to SI units

   intensity=(abs(uxtot).^2+abs(uytot).^2);
   Wio=interp1(ion_lookup(:,1),ion_lookup(:,2),intensity);

   Wio(isnan(Wio))=0;
   net=zeros(length(r),length(t));
   
   for tt=2:length(t)
      net(:,tt)=net(:,tt-1)+Wio(:,tt-1).*(ne0-net(:,tt-1))*dt;
   end
   
%    net=(interp1(tfield,net.',t,'linear','extrap')).';
%    Wio=(interp1(tfield,Wio2.',t,'linear','extrap')).';
   
   absutot2=abs(uxtot.^2)+abs(uytot.^2);
   absutot2loss=absutot2-dz*Wio*ne0*Ip;
   absux=sqrt(abs(absutot2loss.*abs(uxtot.^2)./absutot2));
   absux(isnan(absux))=0;
   absuy=sqrt(abs(absutot2loss.*abs(uytot.^2)./absutot2));
   absuy(isnan(absuy))=0;

   uxtot=absux.*exp(1i*angle(uxtot));
   uytot=absuy.*exp(1i*angle(uytot));
   
   nplasma=-charge.^2/me/eps0/2/w0.^2.*net;
   
   uxtot=uxtot.*exp(-dz*1i*nplasma*w0/c0);
   [Ftx,~]=gradient(1i*nplasma*w0/c0.*uxtot,dt,1);
   uxtot=uxtot-SS_Ion_correction*1i/w0*Ftx*dz;
   
   uytot=uytot.*exp(-dz*1i*nplasma*w0/c0);
   [Fty,~]=gradient(1i*nplasma*w0/c0.*uytot,dt,1);
   uytot=uytot-SS_Ion_correction*1i/w0*Fty*dz;
   
   
   % Projection
   
      for km=1:nmode
        ux(km,:)=squeeze(sum(uxtot.*conj(squeeze(At(km,:,:)))*2*pi.*R*dr,1));
        uy(km,:)=squeeze(sum(uytot.*conj(squeeze(At(km,:,:)))*2*pi.*R*dr,1));      
      end

      
   % Dispersion
   
   uxw=fftshift(fft(ux,[],2),2);
   uxw=uxw.*exp(1i*dz/2.*(-beta+Wm/vg0)).*exp(-alpha*dz/2);
   uxw(:,W<0)=0;
   ux=ifft(ifftshift(uxw,2),[],2);
   
   uyw=fftshift(fft(uy,[],2),2);
   uyw=uyw.*exp(1i*dz/2.*(-beta+Wm/vg0)).*exp(-alpha*dz/2);
   uyw(:,W<0)=0;   
   uy=ifft(ifftshift(uyw,2),[],2);
   
   
   Bint_Kerr=max(max(max(w0/c0*n2*(zobiwanx)*dz)),max(max(w0/c0*n2*(zobiwany)*dz)));%on axis
   Bint_Ion=max(max(abs(nplasma*w0/c0*dz))); %on axis
   %Bint_Ion=0;
   
   
   delta=0.8*min(Bint_Kerr_max/Bint_Kerr,Bint_Ion_max/Bint_Ion);
   
    if Bint_Kerr<Bint_Kerr_max && Bint_Ion<Bint_Ion_max
       index_z=index_z+1;
       z(index_z)=z(index_z-1)+dz;
       uxn=ux;
       uyn=uy;
       ux1(index_z,:)=uxn(1,:);
       uy1(index_z,:)=uyn(1,:);
       
       for km=1:nmode
       Emodex(km,index_z)=sum(sum(abs(ones(nr,1)*uxn(km,:).*squeeze(At(km,:,:))).^2*2*pi.*R*dr*dt));
       Emodey(km,index_z)=sum(sum(abs(ones(nr,1)*uyn(km,:).*squeeze(At(km,:,:))).^2*2*pi.*R*dr*dt));
       end

       Euxtot(index_z)=sum(sum(uxtot.*conj(uxtot)*2*pi.*R*dr*dt));
       Euytot(index_z)=sum(sum(uytot.*conj(uytot)*2*pi.*R*dr*dt));
       Ion_fraction(index_z)=max(net(:,end))/ne0;

       if any(isnan(uxtot))
       disp('DIVERGENCE_uxtot');
       break
       end
       
       Bintv_Kerr(index_z)=Bint_Kerr; %on axis
       Bintv_Ion(index_z)=Bint_Ion;   %on axis

   end
   
   if delta<=0.1
       dz=0.1*dz;
   elseif delta>=4
       dz=4*dz;
   else
       dz=delta*dz;
   end
   
   if dz>dz_max
       dz=dz_max;
   end
   
   if z(index_z)>=Z
       FLAG=1;
   elseif z(index_z)+dz>Z
       dz=Z-z(index_z);
   elseif dz<dz_min
       disp('DIVERGENCE_error');
       break
   end

end

z=z(1:index_z);
ux1=ux1(1:index_z,:);
uy1=uy1(1:index_z,:);

Bintv_Kerr=Bintv_Kerr(1:index_z);
Bintv_Ion=Bintv_Ion(1:index_z);
Euxtot=Euxtot(1:index_z);
Euytot=Euytot(1:index_z);
Emodex=Emodex(:,1:index_z);
Emodey=Emodey(:,1:index_z);
Ion_fraction=Ion_fraction(1:index_z);

uxw1=fftshift(fft(ux1,[],2),2);
uyw1=fftshift(fft(uy1,[],2),2);

uout=ux1(end,:); %lin
%uout=1/sqrt(2)*(ux1(end,:)-1i*uy1(end,:)); %circ
uoutw=fftshift(fft(uout));

E_guided_v=sum(Emodex,1)+sum(Emodey,1);

% Soliton parameters
% Aeff=pi*(0.64*a)^2;
% T0=tgauss/1.76;
% gamma=n2*W(nt/2)/c0/Aeff;
% LD=T0^2/abs(beta2pump);
% LNL=1/gamma/Pcrete;
% Nsol=sqrt(LD/LNL);

% Compression
%
phi2compv=phi2negmax:dphi2:0;
maxi=0;
Ncompopt=1;

for k=1:length(phi2compv)
    
phi2comp=phi2compv(k);
phi3comp=0;
ecompw=fftshift(fft(uout));
ecompw=ecompw.*exp(-1i.*(w.^2/2*phi2comp+w.^3/6*phi3comp));
ecomp=ifft(fftshift(ecompw));
maximoum=max(abs(ecomp.^2));

if maximoum>maxi
Ncompopt=k;
maxi=maximoum;
end

end

phi2compopt=phi2compv(Ncompopt);
phi3compopt=0;
ecompw=fftshift(fft(uout));
ecompw=ecompw.*exp(-1i.*(w.^2/2*phi2compopt+w.^3/6*phi3compopt));
ecomp=ifft(fftshift(ecompw));

% Spectrogram
% 
% windowsize=50e-15;
% windowfunc=window(@hann,round(windowsize/dt));
% windowfunc=[windowfunc.' zeros(1,nt-length(windowfunc))];
% spectrogram=zeros(nt,nt);
% 
%  for zob=1:nt
%      windowshift=circshift(windowfunc,[0 zob-1]);
%      spectrogram(:,zob)=abs(fftshift(fft(uout.*windowshift))).^2;
%  end

% Plots

figure(1);
subplot(211);
plot(z,Emodex.'/energie,z,E_guided_v/energie,'w',z,Euxtot/energie,'w--');
xlabel('Position (m)')
ylabel('Relative energy')

subplot(212);
plot(z,Emodey/energie,z,E_guided_v/energie,'w',z,Euytot/energie,'w--');
xlabel('Position (m)')
ylabel('Relative energy');

figure(2);
subplot(211);
plot(t*1e15,abs(uout).^2./max(abs(uout).^2),t*1e15,abs(ex(1,:).^2)./max(abs(ex(1,:).^2)),t*1e15,abs(ecomp).^2./max(abs(ecomp).^2));
axis([-100 100 0 1.1]);
xlabel('Time (fs)');
ylabel('Intensity (a. u.)');

subplot(212);
plot(lambdav(lambdav>0)*1e9,abs(uxw1(2,lambdav>0)).^2./max(abs(uxw1(2,lambdav>0)).^2),lambdav(lambdav>0)*1e9,abs(uoutw(lambdav>0)).^2./max(abs(uoutw(lambdav>0)).^2));
axis([100 1800 0 1.1])
%plot(w,abs(fftshift(fft(uout)).^2)./max(abs(fftshift(fft(uout)).^2)),w,abs(fftshift(fft(ex1)).^2)./max(abs(fftshift(fft(ex1)).^2)));axis([-1e16 -0.4e16 0 1])
xlabel('Wavelength (nm)')
ylabel('Spectrum (a. u.)')


figure(3);
subplot(211)
plot(z,Bintv_Kerr,z,Bintv_Ion);
legend('Kerr','Ionization');
xlabel('Position (m)');
ylabel('Nonlinear phase (rad)');

subplot(212)
plot(z,Ion_fraction*100);
xlabel('Position (m)');
ylabel('Ionization fraction (%)');
axis([z(1) z(end) 0 max(Ion_fraction*100)*1.1]);

figure(4);
subplot(211);
pcolor(t*1e15,z,abs(ux1.^2)./max(max(abs(ux1.^2))));
axis([-100 100 0 z(end)]);
shading flat
xlabel('Time (fs)')
ylabel('Position (m)')
colorbar;
title('Temporal profile (a. u.)');  
subplot(212);
pcolor(lambdav*1e9,z,10*log10(abs(uxw1.^2)./max(max(abs(uxw1.^2)))));
axis([100 1800 0 z(end)]);
clim([-30 0]);
colorbar;
shading flat
xlabel('Wavelength (nm)')
ylabel('Position (m)')
title('Spectrum (dB)');

figure(5);
pcolor(t*1e15,r*1e6,abs(uxtot.^2));
axis([-100 100 0 a*1e6]);
shading flat;
ylabel('Position (µm)');
xlabel('Time (fs)');

figure(6);
plot(z,gradient(z)*1e3);
xlabel('Position (m)');
ylabel('Stepsize (mm)');

% figure(7)
% plot3(tfield*1e15,efieldx(1,:)./(max(max(abs(efieldx(1,:))),max(abs(efieldy(1,:))))),efieldy(1,:)./(max(max(abs(efieldx(1,:))),max(abs(efieldy(1,:))))),'LineWidth',2);
% hold on
% plot3(tfield*1e15,ones(1,length(tfield)),efieldy(1,:)./(max(max(abs(efieldx(1,:))),max(abs(efieldy(1,:))))));
% plot3(tfield*1e15,efieldx(1,:)./(max(max(abs(efieldx(1,:))),max(abs(efieldy(1,:))))),-1*ones(1,length(tfield)));
% plot3(1000*ones(1,length(tfield)),efieldx(1,:)./(max(max(abs(efieldx(1,:))),max(abs(efieldy(1,:))))),efieldy(1,:)./(max(max(abs(efieldx(1,:))),max(abs(efieldy(1,:))))));
% hold off
% grid on
% xlim([-1000,1000])
% ylim([-1 1]);
% zlim([-1 1]);
% xlabel('Time (fs)')
% ylabel('Ex field (a. u.)')
% zlabel('Ey field (a. u.)');

% figure(6);
% plot(lambdav,10*log10(abs(uxw1(end,:).^2)./max(max(abs(uxw1(end,:).^2)))));
% axis([100e-9 1500e-9 -50 0]);

% figure(7)
% pcolor(t*1e15,lambdav*1e9,spectrogram);shading flat;axis([-50 50 600 1500])
% xlabel('Time (fs)');
% ylabel('Wavelength (nm)');

disp(['Transmission = ' num2str(E_guided_v(end)/energie*100,'%.1f') ' %']);
disp(['Max B-integral per step Kerr= ' num2str(max(Bintv_Kerr),'%.3f') ' rad']);
disp(['Max B-integral per step Ion= ' num2str(max(Bintv_Ion),'%.3f') ' rad']);

% disp(['Ppeak/Pcrit_1bar = ' num2str(Pcrete/Pcrit_1bar,'%.2f') ]);
% disp(['Compression GDD = ' num2str(phi2compopt*1e30,'%.0f') ' fs^2']);

% transmission=E_guided_v(end)/energie;
% zob=[lambdav;uxw1(end,:)./max(abs(uxw1(end,:)))];
% save('Spectrum_40W150kHz_circ.mat','zob','transmission');

% RDW plots
%filtrage UVDW
lambdaDW_start=100e-9;
lambdaDW_stop=400e-9;

uxwDW=uxw1(end,:);
uxwDW(lambdav<lambdaDW_start | lambdav>lambdaDW_stop)=0;
uxDW=ifft(fftshift(uxwDW));
uxDWtot=ones(nr,1)*uxDW.*squeeze(At(1,:,:));
uywDW=uyw1(end,:);
uywDW(lambdav<lambdaDW_start | lambdav>lambdaDW_stop)=0;
uyDW=ifft(fftshift(uywDW));
uyDWtot=ones(nr,1)*uyDW.*squeeze(At(1,:,:));

PDWtot=squeeze(sum((abs(uxDWtot.^2)+abs(uyDWtot.^2))*2*pi.*R*dr,1));
EnergieDWtot=sum(PDWtot*dt);

figure(8)
subplot(211)
plot(t*1e15,abs(uxDW.^2)./max(abs(uxDW.^2)),'r');%,t*1e15,abs(uyDW.^2)./max(abs(uyDW.^2)));
axis([-100 100 0 1]);
xlabel('Time (fs)');
ylabel('UV Intensity (a. u.)')
subplot(212)
%plot(lambdav*1e9,10*log10(abs(uoutw).^2./max(abs(uoutw).^2)),lambdav(lambdav>lambdaDW_start & lambdav<lambdaDW_stop)*1e9,10*log10(abs(uoutw(lambdav>lambdaDW_start & lambdav<lambdaDW_stop)).^2./max(abs(uoutw).^2)),'r');
plot(W/2/pi*1e-12,abs(uoutw).^2./max(abs(uoutw).^2),W(lambdav>lambdaDW_start & lambdav<lambdaDW_stop)/2/pi*1e-12,abs(uoutw(lambdav>lambdaDW_start & lambdav<lambdaDW_stop)).^2./max(abs(uoutw).^2),'r');
%axis([200 1800 -30 0]);
%xlabel('Wavelength (nm)');
%ylabel('Spectrum (dB)');
xlabel('Frequency (THz)');
ylabel('Spectrum (a. u.)');
axis([100 2000 0 1.1]);

% figure(11);
% subplot(221)
% plot(t*1e15,abs(uxDW.^2)./max(abs(uxDW.^2)),'r','LineWidth',2);%t*1e15,abs(uyDW.^2)./max(abs(uyDW.^2)));
% axis([-100 100 0 1]);
% xlabel('Time (fs)');
% ylabel('UV Intensity (a. u.)')
% set(gca,'FontSize',16);
% subplot(222)
% plot(lambdav*1e9,10*log10(abs(uoutw).^2./max(abs(uoutw).^2)),lambdav(lambdav>lambdaDW_start & lambdav<lambdaDW_stop)*1e9,10*log10(abs(uoutw(lambdav>lambdaDW_start & lambdav<lambdaDW_stop)).^2./max(abs(uoutw).^2)),'r','LineWidth',2);
% axis([200 1800 -30 0]);
% xlabel('Wavelength (nm)');
% ylabel('Spectrum (dB)');
% set(gca,'FontSize',16);
% subplot(223);
% pcolor(t*1e15,z,abs(ux1.^2)./max(max(abs(ux1.^2))));
% axis([-100 100 0 z(end)]);
% shading flat
% xlabel('Time (fs)')
% ylabel('Position (m)')
% set(gca,'FontSize',16);
% subplot(224);
% pcolor(lambdav*1e9,z,10*log10(abs(uxw1.^2)./max(max(abs(uxw1.^2)))));
% axis([200 1800 0 z(end)]);
% clim([-30 0]);
% shading flat
% xlabel('Wavelength (nm)')
% ylabel('Position (m)')
% set(gca,'FontSize',16);

% fig = figure(11);
% % setup bottom axis
% ax = axes();
% hold(ax);
% %ax.YAxis.Scale = 'log';
% xlabel(ax, 'Frequency (THz)', 'FontSize', 16);
% ylabel(ax, 'Spectrum (a.u.)', 'FontSize', 16);
% % setup top axis
% ax_top = axes(); % axis to appear at top
% hold(ax_top);
% ax_top.XAxisLocation = 'top';
% ax_top.YAxisLocation = 'right';
% ax_top.YTick = [];
% % ax_top.XDir = 'reverse';
% ax_top.Color = 'none';
% xlabel(ax_top, 'Wavelength (nm)', 'FontSize', 16);
% % linking axis
% linkprop([ax, ax_top],{'Units','Position','ActivePositionProperty'});
% ax.Position(4) = ax.Position(4);
% Freq = [200 400 600 800 1000 1200];%m
% Wavelength = c0./(Freq*1e12)*1e9;
% % configure limits of bottom axis
% ax.XLim = [100 1200];
% ax.XTick = Freq;
% ax.YLim = [0 1.1];
% ax.XAxis.TickLength = [0.02, 0.00];
% ax.YAxis.TickLength = [0.02, 0.00];
% %ax.XAxis.MinorTick = 'off';
% %ax.XAxis.MinorTickValues = linspace(100,1200,19);
% % configure limits and labels of top axis
% y_ticks = [200 300 400 500 600 800 1000 1300 1800];
% Freq_y_tick = c0./(y_ticks*1e-9)*1e-12;
% ax_top.XLim = [Freq(1) Freq(end)];
% ax_top.XTick = fliplr(Freq_y_tick);
% ax_top.XTickLabel = compose('%1.0f', fliplr(y_ticks));
% ax_top.XAxis.TickLength = [0.02, 0.00];
% ax_top.XAxis.MinorTick = 'off';
% ax_top.YLim = [0 1.1];
% 
% %plot(ax, data(:,1), data(:,2), 'r-', 'LineWidth', 3);
% 
% plot(ax,W/2/pi*1e-12,abs(uoutw).^2./max(abs(uoutw).^2),'Linewidth',3);%,W(lambdav>lambdaDW_start & lambdav<lambdaDW_stop)/2/pi*1e-12,abs(uoutw(lambdav>lambdaDW_start & lambdav<lambdaDW_stop)).^2./max(abs(uoutw).^2),'r');
% a=area(ax,W(lambdav>lambdaDW_start & lambdav<lambdaDW_stop)/2/pi*1e-12,abs(uoutw(lambdav>lambdaDW_start & lambdav<lambdaDW_stop)).^2./max(abs(uoutw).^2));
% a.FaceColor = [0.5 0.5 0.5];
% fontsize(fig, 20, "points");
% 
% figure(13)
% pcolor(lambdav*1e9,z,10*log10(abs(uxw1.^2)./max(max(abs(uxw1.^2)))));
% axis([200 1800 0 z(end)]);
% clim([-30 0]);
% shading flat
% xlabel('Wavelength (nm)')
% ylabel('Position (m)')
% set(gca,'FontSize',16);

disp(['Energy in DW = ' num2str(EnergieDWtot*1e6,'%.2f') ' µJ']);
toc;