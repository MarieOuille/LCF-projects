% Propagation in a capillary 2D version for inclusion of nonradially symmetric vector modes
% Instantaneous and delayed (Raman) Kerr, self-steepening, intensity-based ionization
% Vectorial
% Full waveguide and material dispersion
% Pressure gradient, mode resolved in matrix form
% Checked
% Adaptive step based on B integral

% TODO: Choose modes ?

clearvars;
tic;

% Initial condition and capillary geometry
lambda0=1030e-9;
tgauss=30e-15/(sqrt(2*log(2))); %25 fs
energie=120e-6;                 %150e-6

pol='X'; %X or Y or circR or circL
input_type='Gauss'; %Gauss or LG10 or LG20 or Exp_Speck or Exp_Ek
File='spectrum_measured_dscan_yoann.txt';

a=75e-6;
ncap=1.45;
pressurein=4;
pressureout=4;
gas='Ne'; %He or Ne or Ar or Xe or N2 !Warning! tables ionisations pas toutes dispo pour le moment 
Z=1;

% Temporal and spatial ranges and steps

dt=0.2e-15;
T=500e-15;

dz_min=1e-6;
dz_max=50e-3;
Bint_Kerr_max=0.1;
Bint_Ion_max=0.1;
nz_max=10000;

dx=5e-6;

% Compression parameters
phi2negmax=-20000e-30;
dphi2=100e-30;

%_____________________________________________________

% Parameters
c0=299792458;
charge=1.602e-19;
me=9.11e-31;
hbar=6.62e-34;
eps0=8.85e-12;
mu0=4*pi*1e-7;
w0=2*pi*c0/lambda0;
rgauss=0.64*a;
rOAM1=0.56*a;
rOAM2=0.52*a;

% Time/frequency variables
%dt=1.01*lambda0/2/c0;
nt=2^(ceil(log2(T/dt)));
t=-(nt/2-1)*dt:dt:(nt/2)*dt;
w=2*pi*(1/dt)*(-nt/2:nt/2-1)./nt;
dw=w(2)-w(1);

lambdav=c0*2*pi./(w+w0);
lmicron=lambdav*1e6;
W=w+w0;

% Transverse spatial variable
dx=a/round(a/dx);
x=-a:dx:a;
nx=length(x);
y=x;

% Capillary modes

[X,Y,T]=ndgrid(x,y,t);
R=sqrt(X.^2+Y.^2);
TETA=atan2(Y,X);

[X2,Y2]=ndgrid(x,y); 
R2=sqrt(X2.^2+Y2.^2);
TETA2=atan2(Y2,X2);

u01=2.40482555769577;
u02=5.52007811028631;
u03=8.65372791291101;
u11=3.83170597020751;
u12=7.01558666981562;
u21=5.135622301840683;

nmode=18;
Ax=zeros(nmode,nx,nx);
Axt=ones(nmode,nx,nx,nt);
Ay=zeros(nmode,nx,nx);
Ayt=zeros(nmode,nx,nx,nt);
alpha=zeros(nmode,nt);
Wm=ones(nmode,1)*w;
ukm=zeros(1,nmode);

%mode 1: EH11x
ukm(1)=u01;
Ax(1,:,:)=besselj(0,ukm(1)*R2/a);
Ax(1,R2>a)=0;
Ay(1,:,:)=0;
norm=sqrt(sum(sum(abs(Ax(1,:,:).^2).*dx*dx))+sum(sum(abs(Ay(1,:,:).^2).*dx*dx)));

Ax(1,:,:)=Ax(1,:,:)./norm;
Ay(1,:,:)=Ay(1,:,:)./norm; 

Axt(1,:,:,:)=repmat(squeeze(Ax(1,:,:)),1,1,nt);
Ayt(1,:,:,:)=repmat(squeeze(Ay(1,:,:)),1,1,nt);

alpha(1,:)=(ukm(1)/2/pi)^2*lambdav.^2/a^3*(0.5*(ncap^2+1)/sqrt(ncap^2-1));

%mode 2: EH11y
ukm(2)=u01;
Ax(2,:,:)=0;
Ay(2,:,:)=besselj(0,ukm(2)*R2/a);
Ay(2,R2>a)=0;
norm=sqrt(sum(sum(abs(Ax(2,:,:).^2).*dx*dx))+sum(sum(abs(Ay(2,:,:).^2).*dx*dx)));

Ax(2,:,:)=Ax(2,:,:)./norm;
Ay(2,:,:)=Ay(2,:,:)./norm; 

Ayt(2,:,:,:)=repmat(squeeze(Ay(2,:,:)),1,1,nt);
Axt(2,:,:,:)=repmat(squeeze(Ax(2,:,:)),1,1,nt);

alpha(2,:)=(ukm(2)/2/pi)^2*lambdav.^2/a^3*(0.5*(ncap^2+1)/sqrt(ncap^2-1));

%mode 3: EH12x

ukm(3)=u02;
Ax(3,:,:)=besselj(0,ukm(3)*R2/a);
Ax(3,R2>a)=0;
Ay(3,:,:)=0;
norm=sqrt(sum(sum(abs(Ax(3,:,:).^2).*dx*dx))+sum(sum(abs(Ay(3,:,:).^2).*dx*dx)));

Ax(3,:,:)=Ax(3,:,:)./norm;
Ay(3,:,:)=Ay(3,:,:)./norm; 

Axt(3,:,:,:)=repmat(squeeze(Ax(3,:,:)),1,1,nt);
Ayt(3,:,:,:)=repmat(squeeze(Ay(3,:,:)),1,1,nt);

alpha(3,:)=(ukm(3)/2/pi)^2*lambdav.^2/a^3*(0.5*(ncap^2+1)/sqrt(ncap^2-1));

%mode 4: EH12y
ukm(4)=u02;
Ax(4,:,:)=0;
Ay(4,:,:)=besselj(0,ukm(4)*R2/a);
Ay(4,R2>a)=0;
norm=sqrt(sum(sum(abs(Ax(4,:,:).^2).*dx*dx))+sum(sum(abs(Ay(4,:,:).^2).*dx*dx)));

Ax(4,:,:)=Ax(4,:,:)./norm;
Ay(4,:,:)=Ay(4,:,:)./norm; 

Ayt(4,:,:,:)=repmat(squeeze(Ay(4,:,:)),1,1,nt);
Axt(4,:,:,:)=repmat(squeeze(Ax(4,:,:)),1,1,nt);

alpha(4,:)=(ukm(4)/2/pi)^2*lambdav.^2/a^3*(0.5*(ncap^2+1)/sqrt(ncap^2-1));

%mode 5: EH21x
ukm(5)=u11;
Ateta=besselj(1,ukm(5)*R2/a).*cos(2*TETA2);
Ar=besselj(1,ukm(5)*R2/a).*sin(2*TETA2);

Ax(5,:,:)=Ateta.*-sin(TETA2)+Ar.*cos(TETA2);
Ax(5,R2>a)=0;
Ay(5,:,:)=Ateta.*cos(TETA2)+Ar.*sin(TETA2);
Ay(5,R2>a)=0;

norm=sqrt(sum(sum(abs(Ax(5,:,:).^2).*dx*dx))+sum(sum(abs(Ay(5,:,:).^2).*dx*dx)));

Ax(5,:,:)=Ax(5,:,:)./norm;
Ay(5,:,:)=Ay(5,:,:)./norm; 

Axt(5,:,:,:)=repmat(squeeze(Ax(5,:,:)),1,1,nt);
Ayt(5,:,:,:)=repmat(squeeze(Ay(5,:,:)),1,1,nt);

alpha(5,:)=(ukm(5)/2/pi)^2*lambdav.^2/a^3*(0.5*(ncap^2+1)/sqrt(ncap^2-1));

%mode 6: EH21y
ukm(6)=u11;
Ateta=besselj(1,ukm(6)*R2/a).*cos(2*TETA2-pi/2);
Ar=besselj(1,ukm(6)*R2/a).*sin(2*TETA2-pi/2);

Ax(6,:,:)=Ateta.*-sin(TETA2)+Ar.*cos(TETA2);
Ax(6,R2>a)=0;
Ay(6,:,:)=Ateta.*cos(TETA2)+Ar.*sin(TETA2);
Ay(6,R2>a)=0;

norm=sqrt(sum(sum(abs(Ax(6,:,:).^2).*dx*dx))+sum(sum(abs(Ay(6,:,:).^2).*dx*dx)));

Ax(6,:,:)=Ax(6,:,:)./norm;
Ay(6,:,:)=Ay(6,:,:)./norm; 

Axt(6,:,:,:)=repmat(squeeze(Ax(6,:,:)),1,1,nt);
Ayt(6,:,:,:)=repmat(squeeze(Ay(6,:,:)),1,1,nt);

alpha(6,:)=(ukm(6)/2/pi)^2*lambdav.^2/a^3*(0.5*(ncap^2+1)/sqrt(ncap^2-1));

%mode 7 : TE01
ukm(7)=u11;
Ateta=besselj(1,ukm(7)*R2/a);

Ax(7,:,:)=Ateta.*-sin(TETA2);
Ax(7,R2>a)=0;
Ay(7,:,:)=Ateta.*cos(TETA2);
Ay(7,R2>a)=0;

norm=sqrt(sum(sum(abs(Ax(7,:,:).^2).*dx*dx))+sum(sum(abs(Ay(7,:,:).^2).*dx*dx)));

Ax(7,:,:)=Ax(7,:,:)./norm;
Ay(7,:,:)=Ay(7,:,:)./norm; 

Axt(7,:,:,:)=repmat(squeeze(Ax(7,:,:)),1,1,nt);
Ayt(7,:,:,:)=repmat(squeeze(Ay(7,:,:)),1,1,nt);

alpha(7,:)=(ukm(7)/2/pi)^2*lambdav.^2/a^3*(1/sqrt(ncap^2-1));

%mode 8 : TM01
ukm(8)=u11;
Ar=besselj(1,ukm(8)*R2/a);

Ax(8,:,:)=Ar.*cos(TETA2);
Ax(8,R2>a)=0;
Ay(8,:,:)=Ar.*sin(TETA2);
Ay(8,R2>a)=0;

norm=sqrt(sum(sum(abs(Ax(8,:,:).^2).*dx*dx))+sum(sum(abs(Ay(8,:,:).^2).*dx*dx)));

Ax(8,:,:)=Ax(8,:,:)./norm;
Ay(8,:,:)=Ay(8,:,:)./norm; 

Axt(8,:,:,:)=repmat(squeeze(Ax(8,:,:)),1,1,nt);
Ayt(8,:,:,:)=repmat(squeeze(Ay(8,:,:)),1,1,nt);

alpha(8,:)=(ukm(8)/2/pi)^2*lambdav.^2/a^3*(ncap^2/sqrt(ncap^2-1));

%mode 9 : EH22x
ukm(9)=u12;
Ateta=besselj(1,ukm(9)*R2/a).*cos(2*TETA2);
Ar=besselj(1,ukm(9)*R2/a).*sin(2*TETA2);

Ax(9,:,:)=Ateta.*-sin(TETA2)+Ar.*cos(TETA2);
Ax(9,R2>a)=0;
Ay(9,:,:)=Ateta.*cos(TETA2)+Ar.*sin(TETA2);
Ay(9,R2>a)=0;

norm=sqrt(sum(sum(abs(Ax(9,:,:).^2).*dx*dx))+sum(sum(abs(Ay(9,:,:).^2).*dx*dx)));

Ax(9,:,:)=Ax(9,:,:)./norm;
Ay(9,:,:)=Ay(9,:,:)./norm; 

Axt(9,:,:,:)=repmat(squeeze(Ax(9,:,:)),1,1,nt);
Ayt(9,:,:,:)=repmat(squeeze(Ay(9,:,:)),1,1,nt);

alpha(9,:)=(ukm(9)/2/pi)^2*lambdav.^2/a^3*(0.5*(ncap^2+1)/sqrt(ncap^2-1));

%mode 10: EH22y
ukm(10)=u12;
Ateta=besselj(1,ukm(10)*R2/a).*cos(2*TETA2-pi/2);
Ar=besselj(1,ukm(10)*R2/a).*sin(2*TETA2-pi/2);

Ax(10,:,:)=Ateta.*-sin(TETA2)+Ar.*cos(TETA2);
Ax(10,R2>a)=0;
Ay(10,:,:)=Ateta.*cos(TETA2)+Ar.*sin(TETA2);
Ay(10,R2>a)=0;

norm=sqrt(sum(sum(abs(Ax(10,:,:).^2).*dx*dx))+sum(sum(abs(Ay(10,:,:).^2).*dx*dx)));

Ax(10,:,:)=Ax(10,:,:)./norm;
Ay(10,:,:)=Ay(10,:,:)./norm; 

Axt(10,:,:,:)=repmat(squeeze(Ax(10,:,:)),1,1,nt);
Ayt(10,:,:,:)=repmat(squeeze(Ay(10,:,:)),1,1,nt);

alpha(10,:)=(ukm(10)/2/pi)^2*lambdav.^2/a^3*(0.5*(ncap^2+1)/sqrt(ncap^2-1));

%mode 11 : TE02
ukm(11)=u12;
Ateta=besselj(1,ukm(11)*R2/a);

Ax(11,:,:)=Ateta.*-sin(TETA2);
Ax(11,R2>a)=0;
Ay(11,:,:)=Ateta.*cos(TETA2);
Ay(11,R2>a)=0;

norm=sqrt(sum(sum(abs(Ax(11,:,:).^2).*dx*dx))+sum(sum(abs(Ay(11,:,:).^2).*dx*dx)));

Ax(11,:,:)=Ax(11,:,:)./norm;
Ay(11,:,:)=Ay(11,:,:)./norm; 

Axt(11,:,:,:)=repmat(squeeze(Ax(11,:,:)),1,1,nt);
Ayt(11,:,:,:)=repmat(squeeze(Ay(11,:,:)),1,1,nt);

alpha(11,:)=(ukm(11)/2/pi)^2*lambdav.^2/a^3*(1/sqrt(ncap^2-1));

%mode 12 : TM02
ukm(12)=u12;
Ar=besselj(1,ukm(12)*R2/a);

Ax(12,:,:)=Ar.*cos(TETA2);
Ax(12,R2>a)=0;
Ay(12,:,:)=Ar.*sin(TETA2);
Ay(12,R2>a)=0;

norm=sqrt(sum(sum(abs(Ax(12,:,:).^2).*dx*dx))+sum(sum(abs(Ay(12,:,:).^2).*dx*dx)));

Ax(12,:,:)=Ax(12,:,:)./norm;
Ay(12,:,:)=Ay(12,:,:)./norm; 

Axt(12,:,:,:)=repmat(squeeze(Ax(12,:,:)),1,1,nt);
Ayt(12,:,:,:)=repmat(squeeze(Ay(12,:,:)),1,1,nt);

alpha(12,:)=(ukm(12)/2/pi)^2*lambdav.^2/a^3*(ncap^2/sqrt(ncap^2-1));

%mode 13: EH13x

ukm(13)=u03;
Ax(13,:,:)=besselj(0,ukm(13)*R2/a);
Ax(13,R2>a)=0;
Ay(13,:,:)=0;
norm=sqrt(sum(sum(abs(Ax(13,:,:).^2).*dx*dx))+sum(sum(abs(Ay(13,:,:).^2).*dx*dx)));

Ax(13,:,:)=Ax(13,:,:)./norm;
Ay(13,:,:)=Ay(13,:,:)./norm; 

Axt(13,:,:,:)=repmat(squeeze(Ax(13,:,:)),1,1,nt);
Ayt(13,:,:,:)=repmat(squeeze(Ay(13,:,:)),1,1,nt);

alpha(13,:)=(ukm(13)/2/pi)^2*lambdav.^2/a^3*(0.5*(ncap^2+1)/sqrt(ncap^2-1));

%mode 14: EH13y
ukm(14)=u03;
Ax(14,:,:)=0;
Ay(14,:,:)=besselj(0,ukm(14)*R2/a);
Ay(14,R2>a)=0;
norm=sqrt(sum(sum(abs(Ax(14,:,:).^2).*dx*dx))+sum(sum(abs(Ay(14,:,:).^2).*dx*dx)));

Ax(14,:,:)=Ax(14,:,:)./norm;
Ay(14,:,:)=Ay(14,:,:)./norm; 

Ayt(14,:,:,:)=repmat(squeeze(Ay(14,:,:)),1,1,nt);
Axt(14,:,:,:)=repmat(squeeze(Ax(14,:,:)),1,1,nt);

alpha(14,:)=(ukm(14)/2/pi)^2*lambdav.^2/a^3*(0.5*(ncap^2+1)/sqrt(ncap^2-1));

%mode 15: EH31x
ukm(15)=u21;
Ateta=besselj(2,ukm(15)*R2/a).*cos(3*TETA2);
Ar=besselj(2,ukm(15)*R2/a).*sin(3*TETA2);

Ax(15,:,:)=Ateta.*-sin(TETA2)+Ar.*cos(TETA2);
Ax(15,R2>a)=0;
Ay(15,:,:)=Ateta.*cos(TETA2)+Ar.*sin(TETA2);
Ay(15,R2>a)=0;

norm=sqrt(sum(sum(abs(Ax(15,:,:).^2).*dx*dx))+sum(sum(abs(Ay(15,:,:).^2).*dx*dx)));

Ax(15,:,:)=Ax(15,:,:)./norm;
Ay(15,:,:)=Ay(15,:,:)./norm; 

Axt(15,:,:,:)=repmat(squeeze(Ax(15,:,:)),1,1,nt);
Ayt(15,:,:,:)=repmat(squeeze(Ay(15,:,:)),1,1,nt);

alpha(15,:)=(ukm(15)/2/pi)^2*lambdav.^2/a^3*(0.5*(ncap^2+1)/sqrt(ncap^2-1));

%mode 16: EH31y
ukm(16)=u21;
Ateta=besselj(2,ukm(16)*R2/a).*cos(3*TETA2-pi/2);
Ar=besselj(2,ukm(16)*R2/a).*sin(3*TETA2-pi/2);

Ax(16,:,:)=Ateta.*-sin(TETA2)+Ar.*cos(TETA2);
Ax(16,R2>a)=0;
Ay(16,:,:)=Ateta.*cos(TETA2)+Ar.*sin(TETA2);
Ay(16,R2>a)=0;

norm=sqrt(sum(sum(abs(Ax(16,:,:).^2).*dx*dx))+sum(sum(abs(Ay(16,:,:).^2).*dx*dx)));

Ax(16,:,:)=Ax(16,:,:)./norm;
Ay(16,:,:)=Ay(16,:,:)./norm; 

Axt(16,:,:,:)=repmat(squeeze(Ax(16,:,:)),1,1,nt);
Ayt(16,:,:,:)=repmat(squeeze(Ay(16,:,:)),1,1,nt);

alpha(16,:)=(ukm(16)/2/pi)^2*lambdav.^2/a^3*(0.5*(ncap^2+1)/sqrt(ncap^2-1));

%mode 17: EH-11x
ukm(17)=u21;
Ateta=besselj(-2,ukm(17)*R2/a).*cos(-TETA2);
Ar=besselj(-2,ukm(17)*R2/a).*sin(-TETA2);

Ax(17,:,:)=Ateta.*-sin(TETA2)+Ar.*cos(TETA2);
Ax(17,R2>a)=0;
Ay(17,:,:)=Ateta.*cos(TETA2)+Ar.*sin(TETA2);
Ay(17,R2>a)=0;

norm=sqrt(sum(sum(abs(Ax(17,:,:).^2).*dx*dx))+sum(sum(abs(Ay(17,:,:).^2).*dx*dx)));

Ax(17,:,:)=Ax(17,:,:)./norm;
Ay(17,:,:)=Ay(17,:,:)./norm; 

Axt(17,:,:,:)=repmat(squeeze(Ax(17,:,:)),1,1,nt);
Ayt(17,:,:,:)=repmat(squeeze(Ay(17,:,:)),1,1,nt);

alpha(17,:)=(ukm(17)/2/pi)^2*lambdav.^2/a^3*(0.5*(ncap^2+1)/sqrt(ncap^2-1));

%mode 18: EH-11y
ukm(18)=u21;
Ateta=besselj(-2,ukm(18)*R2/a).*cos(-TETA2-pi/2);
Ar=besselj(-2,ukm(18)*R2/a).*sin(-TETA2-pi/2);

Ax(18,:,:)=Ateta.*-sin(TETA2)+Ar.*cos(TETA2);
Ax(18,R2>a)=0;
Ay(18,:,:)=Ateta.*cos(TETA2)+Ar.*sin(TETA2);
Ay(18,R2>a)=0;

norm=sqrt(sum(sum(abs(Ax(18,:,:).^2).*dx*dx))+sum(sum(abs(Ay(18,:,:).^2).*dx*dx)));

Ax(18,:,:)=Ax(18,:,:)./norm;
Ay(18,:,:)=Ay(18,:,:)./norm; 

Axt(18,:,:,:)=repmat(squeeze(Ax(18,:,:)),1,1,nt);
Ayt(18,:,:,:)=repmat(squeeze(Ay(18,:,:)),1,1,nt);

alpha(18,:)=(ukm(18)/2/pi)^2*lambdav.^2/a^3*(0.5*(ncap^2+1)/sqrt(ncap^2-1));




%______________________


alpha=real(alpha);

% Input field

switch input_type
    case 'LG10' 
    Ex=R/rOAM1.*exp(1i*TETA).*exp(-R.^2/rOAM1.^2).*exp(-(T/tgauss).^2);
    energieint=sum(sum(sum(abs(Ex.^2))))*dt*dx*dx;
    Ex=sqrt(energie/energieint).*Ex;
    
    case 'LG20'
    Ex=(R/rOAM2).^2.*exp(1i*2*TETA).*exp(-R.^2/rOAM2.^2).*exp(-(T/tgauss).^2);
    energieint=sum(sum(sum(abs(Ex.^2))))*dt*dx*dx;
    Ex=sqrt(energie/energieint).*Ex;
   
    case 'Gauss'
    Ex=exp(-(R/rgauss).^2).*exp(-(T/tgauss).^2);
    energieint=sum(sum(sum(abs(Ex.^2))))*dt*dx*dx;
    Ex=sqrt(energie/energieint).*Ex;

    case 'Exp_Ek'
    aa=load(File);
    tdata=aa(:,1)*1e-15;
    edata=sqrt(aa(:,2)).*exp(1i*aa(:,3));
    einterp=interp1(tdata,edata,t);
    einterp(isnan(einterp))=0;
    
    Ex=ones(nx,nx)*einterp.*exp(-(R/rgauss).^2);
    energieint=sum(sum(sum(abs(Ex.^2))))*dt*dx*dx;
    Ex=sqrt(energie/energieint).*Ex;
    
    case 'Exp_Speck'
    energie=energie*(1-exp(-2*a^2/rgauss^2)); % No energy outside capillary
    aaa=load(File);
    ldata=aaa(:,1)*1e-9;
    aaa(:,2) = aaa(:,2)./max(aaa(:,2));
    ewdata=sqrt(aaa(:,2)).*exp(-1i*aaa(:,3)); %might need to change phase sign
    ewdata=circshift(ewdata,0);
    ewinterp=interp1(ldata,ewdata,lambdav);
    ewinterp(isnan(ewinterp))=0;
    einterp=fftshift(ifft(ifftshift(ewinterp)));
    Ex=repmat(einterp,1,nx*nx);
    Ex=reshape(Ex,nt,nx,nx);
    Ex=permute(Ex,[3 2 1]);
    Ex=Ex.*exp(-(R/rgauss).^2);
    energieint=sum(sum(sum(abs(Ex.^2))))*dt*dx*dx;
    Ex=sqrt(energie/energieint).*Ex;
end

Pcrete=max(squeeze(sum(sum(abs(Ex.^2)*dx*dx,1))));
Icrete=max(max(max(abs(Ex.^2))));

% Polarisation state X and Y components
switch pol
    case 'X'
        ell=1e4;
    case 'Y'
        ell=1e-4;
    case 'circL'
        ell=1;
    case 'circR'
        ell=-1;
end

Ex=ell/sqrt(1+ell^2)*Ex;
Ey=1i/ell*Ex;

e=zeros(nmode,nt);

for km=1:nmode
e(km,:)=squeeze(sum(sum((Ex.*conj(squeeze(Axt(km,:,:,:)))+Ey.*conj(squeeze(Ayt(km,:,:,:))))*dx*dx,1),2));
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
    if or(strcmp(pol,'X'),strcmp(pol,'Y'))
        ion_lookup=load('linear_rate_Ar_m=0_lambda=1030nm.dat');
    elseif or(strcmp(pol,'circL'),strcmp(pol,'circR'))
        ion_lookup=load('circular_rate_Ar_m=0_lambda=1030nm.dat');
    end
    
    
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
    if or(strcmp(pol,'X'),strcmp(pol,'Y'))
        ion_lookup=load('linear_rate_Ne_m=0_lambda=1030nm.dat');
    elseif or(strcmp(pol,'circL'),strcmp(pol,'circR'))
        ion_lookup=load('circular_rate_Ne_m=0.dat');
    end
    
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
    if or(strcmp(pol,'X'),strcmp(pol,'Y'))
        ion_lookup=load('linear_rate_Xe_m=0.dat');
    elseif or(strcmp(pol,'circL'),strcmp(pol,'circR'))
        ion_lookup=load('circular_rate_Xe_m=0.dat');
    end

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
    if or(strcmp(pol,'X'),strcmp(pol,'Y'))
        ion_lookup=load('linear_rate_He_3_lambda=1030nm.dat');
    elseif or(strcmp(pol,'circL'),strcmp(pol,'circR'))
        ion_lookup=load('circular_rate_He_2.dat');
    end

    case 'N2'
    nlambda=(6.8552e-5+3.243157e-2./(144-lmicron.^-2)); %N2
    nlambda(isnan(nlambda))=0;
    nlambda=real(nlambda);
    n2_1bar=4e-23; %N2 global
    tau1=66e-15; %62.5e-15 N2
    tau2=133e-15;  %120e-15 N2        
    fr=0.68;       %fr_1=0.75 for N2
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
end

Pcrit_1bar=lambdav(nt/2)^2/8/n2_1bar;
hr=(tau1^2+tau2^2)/(tau1*tau2^2).*exp(-t/tau2).*sin(t/tau1);
hr(t<0)=0;
hr=hr/sum(hr);
Resp=fr*hr;
Resp(t==0)=Resp(t==0)+(1-fr);
Respw=fft(Resp);
zob=zeros(1,1,length(t));
zob(1,1,:)=Respw;
Respw=repmat(zob,nx,nx);


% Propagation

un=e;

u1=zeros(nz_max,nt);

Emode=zeros(nmode,nz_max);
Emode_test=zeros(nmode,nz_max);

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
    beta(km,:)=W./c0.*sqrt(complex(n.^2-ukm(km)^2*c0^2/a^2./W.^2));
end

beta=real(beta);
       
switch input_type
       case 'Exp_Speck'
       beta1=gradient(beta(1,:),W);
       case 'Gauss'
       beta1=gradient(beta(1,:),W);
       case 'LG10'
       beta1=gradient(beta(5,:),W);
       case 'LG20'
       beta1=gradient(beta(15,:),W);
end

vg0=1/beta1(nt/2);

% Propagation over dz

   uw=fftshift(fft(un,[],2),2);
   uw=uw.*Spectral_window;
   
   % Dispersion
   
   uw=uw.*exp(1i*dz/2.*(-beta+Wm/vg0)).*exp(-alpha*dz/2);
   uw(:,W<0)=0;
   u=ifft(ifftshift(uw,2),[],2);
      
   u=u.*Temporal_window;
   
   % Total field
  
   uxtot=zeros(nx,nx,nt);
   uytot=zeros(nx,nx,nt);
   
   for km=1:nmode

       zob(1,1,:)=squeeze(u(km,:));
       zob3=repmat(zob,nx,nx);
       
       uxtot=uxtot+zob3.*squeeze(Axt(km,:,:,:));
       uytot=uytot+zob3.*squeeze(Ayt(km,:,:,:));
   end
   
   
   % Kerr + Raman nonlinearity
   
   zobiwanx=circshift(fftshift(ifft(fft(abs(uxtot.^2)+2/3*abs(uytot.^2),[],3).*Respw,[],3),3),1,3); %Kerr + Raman
   zobiwany=circshift(fftshift(ifft(fft(abs(uytot.^2)+2/3*abs(uxtot.^2),[],3).*Respw,[],3),3),1,3); %Kerr + Raman
   
   uxtot=uxtot.*exp(-dz*1i*w0/c0*n2*zobiwanx);
   uytot=uytot.*exp(-dz*1i*w0/c0*n2*zobiwany);
   
   [~,~,Ftx]=gradient(+1i*w0/c0*n2*(zobiwanx.*uxtot+1/3*uytot.^2.*conj(uxtot)),dt);
   Ftx(isnan(Ftx))=0;
   uxtot=uxtot+1i/w0*dz*Ftx-1i/3.*w0/c0*n2*uytot.^2.*conj(uxtot)*dz;
   
   [~,~,Fty]=gradient(+1i*w0/c0*n2*(zobiwany.*uytot+1/3*uxtot.^2.*conj(uytot)),dt);
   Fty(isnan(Fty))=0;
   uytot=uytot+1i/w0*dz*Fty-1i/3.*w0/c0*n2*uxtot.^2.*conj(uytot)*dz; 
   
%  Ionization

   %absefield=sqrt((abs(uxtot).^2+abs(uytot).^2)*2/eps0/c0)*2/pi; % cycle-averaged magnitude of the field in V/m (linear polarization)  
                                                                 % underestimates ionization?
                                                                 
   %Wio=kappaio*(muio./absefield).^(2*nstar-1).*exp(-muio/3./absefield); % ADK these malvache
   
%    F_TL= absefield/5.14220652e11; % field strength (atomic units);
%    Wio=C_TL^2 /2^abs(m_TL) /factorial(abs(m_TL))...% Tong Lin 2005 d'après HHGmax
%      * (2*l_TL+1)*factorial(l_TL+abs(m_TL)) /2 /factorial(l_TL-abs(m_TL))...
%      * 1 / kappa_TL^( 2*Z_TL/kappa_TL - 1 )...
%      * (2*kappa_TL^3 ./ F_TL) .^ ( 2*Z_TL/kappa_TL - abs(m_TL) - 1 )...
%     .* exp( -2/3 * kappa_TL^3./F_TL )...
%     .* exp( -alpha_TL * (Z_TL^2/Ip_TL) * (F_TL/kappa_TL^3) );
%    Wio=Wio*(1/2.418884326505e-17); %back to SI units

   intensity=(abs(uxtot).^2+abs(uytot).^2);
   Wio=interp1(ion_lookup(:,1),ion_lookup(:,2),intensity);

   Wio(isnan(Wio))=0;
   net=zeros(nx,nx,nt);
   
   for tt=2:length(t)
      net(:,:,tt)=net(:,:,tt-1)+Wio(:,:,tt-1).*(ne0-net(:,:,tt-1))*dt;
   end
   
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
   [~,~,Ftx]=gradient(1i*nplasma*w0/c0.*uxtot,dt);
   Ftx(isnan(Ftx))=0;
   uxtot=uxtot-1i/w0*Ftx*dz;
   
   uytot=uytot.*exp(-dz*1i*nplasma*w0/c0);
   [~,~,Fty]=gradient(1i*nplasma*w0/c0.*uytot,dt);
   Fty(isnan(Fty))=0;
   uytot=uytot-1i/w0*Fty*dz;
   
   % Projection
   
      for km=1:nmode
        u(km,:)=squeeze(sum(sum((uxtot.*conj(squeeze(Axt(km,:,:,:)))+uytot.*conj(squeeze(Ayt(km,:,:,:))))*dx*dx,1),2));
      end

   % Dispersion
      
   uw=fftshift(fft(u,[],2),2);
   uw=uw.*exp(1i*dz/2.*(-beta+Wm/vg0)).*exp(-alpha*dz/2);
   uw(:,W<0)=0;
   u=ifft(ifftshift(uw,2),[],2);
   
   Bint_Kerr=max(max(max(max(w0/c0*n2*(zobiwanx)*dz))),max(max(max(w0/c0*n2*(zobiwany)*dz))));%on axis
   Bint_Ion=max(max(max(abs(nplasma*w0/c0*dz)))); %on axis
   %Bint_Ion=0;
   
   
   delta=0.8*min(Bint_Kerr_max/Bint_Kerr,Bint_Ion_max/Bint_Ion);
   
    if Bint_Kerr<Bint_Kerr_max && Bint_Ion<Bint_Ion_max
       index_z=index_z+1;
       z(index_z)=z(index_z-1)+dz;
       un=u;
       switch input_type
           case 'Exp_Speck'
                switch pol
                    case 'X'
                        u1(index_z,:)=un(1,:);
                    case 'Y'
                        u1(index_z,:)=un(2,:);
                    case 'circL'
                        u1(index_z,:)=1/sqrt(2)*(un(1,:)-1i*un(2,:));
                    case 'circR'
                        u1(index_z,:)=1/sqrt(2)*(un(1,:)+1i*un(2,:));
                end
           case 'Gauss'
                switch pol
                    case 'X'
                        u1(index_z,:)=un(1,:);
                    case 'Y'
                        u1(index_z,:)=un(2,:);
                    case 'circL'
                        u1(index_z,:)=1/sqrt(2)*(un(1,:)-1i*un(2,:));
                    case 'circR'
                        u1(index_z,:)=1/sqrt(2)*(un(1,:)+1i*un(2,:));
                end
           case 'LG10'
                switch pol
                    case 'X'
                        u1(index_z,:)=un(5,:);  
                    case 'Y'
                        u1(index_z,:)=un(5,:);
                    case 'circL'
                        u1(index_z,:)=un(5,:);
                    case 'circR'
                        u1(index_z,:)=un(7,:);
                end
                
                case 'LG20'
                switch pol
                    case 'X'
                        u1(index_z,:)=un(15,:);  
                    case 'Y'
                        u1(index_z,:)=un(15,:);
                    case 'circL'
                        u1(index_z,:)=un(15,:);
                    case 'circR'
                        u1(index_z,:)=un(17,:);
                end
       end
       
       for km=1:nmode           
       Emode(km,index_z)=sum(abs(u(km,:).^2)*dt);
       end

       Euxtot(index_z)=sum(sum(sum(uxtot.*conj(uxtot)*dx*dx*dt)));
       Euytot(index_z)=sum(sum(sum(uytot.*conj(uytot)*dx*dx*dt)));
       Ion_fraction(index_z)=max(max(net(:,:,end)))/ne0;
       
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
u1=u1(1:index_z,:);

Bintv_Kerr=Bintv_Kerr(1:index_z);
Bintv_Ion=Bintv_Ion(1:index_z);
Euxtot=Euxtot(1:index_z);
Euytot=Euytot(1:index_z);
Emode=Emode(:,1:index_z);
Ion_fraction=Ion_fraction(1:index_z);

uw1=fftshift(fft(u1,[],2),2);
uout=u1(end,:);
uoutw=fftshift(fft(uout));

E_guided_v=sum(Emode,1);

% Soliton parameters
% Aeff=pi*(0.64*a)^2;
% T0=tgauss/1.76;
% gamma=n2*W(nt/2)/c0/Aeff;
% LD=T0^2/abs(beta2pump);
% LNL=1/gamma/Pcrete;
% Nsol=sqrt(LD/LNL);

% Compression

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
tiledlayout(5,1);
nexttile([2 1]);
plot(z,Emode.'/energie,z,E_guided_v/energie,'k',z,Euxtot/energie,'r--',z,Euytot/energie,'b--');
xlabel('Position (m)');
ylabel('Relative energy');
nexttile
plot(z,Bintv_Kerr,z,Bintv_Ion);
legend('Kerr','Ionization');
xlabel('Position (m)');
ylabel('Nonlinear phase per step (rad)');
nexttile
plot(z,gradient(z)*1e3);
xlabel('Position (m)');
ylabel('Stepsize (mm)');
nexttile
plot(z,Ion_fraction*100);
xlabel('Position (m)');
ylabel('Ionization fraction (%)');
axis([z(1) z(end) 0 max(Ion_fraction*100)*1.1]);


figure(2);
subplot(211);
%plot(t*1e15,abs(u1(2,:).^2)./max(abs(u1(2,:).^2)),t*1e15,squeeze(sum(sum(abs(uytot.^2),1),2))./max(squeeze(sum(sum(abs(uytot.^2),1),2))))%,t*1e15,abs(ecomp).^2./max(abs(ecomp).^2));
plot(t*1e15,abs(u1(2,:).^2)./max(abs(u1(2,:).^2)),t*1e15,abs(u1(end,:).^2)./max(abs(u1(end,:).^2)),t*1e15,abs(ecomp).^2./max(abs(ecomp).^2));
axis([t(1)*1e15 t(end)*1e15 0 1.1]);
xlabel('Time (fs)');
ylabel('Intensity (a. u.)');
subplot(212);
plot(lambdav(lambdav>0)*1e9,abs(uw1(2,lambdav>0)).^2./max(abs(uw1(2,lambdav>0)).^2),lambdav(lambdav>0)*1e9,abs(uoutw(lambdav>0)).^2./max(abs(uoutw(lambdav>0)).^2));
axis([100 1800 0 1.1])
xlabel('Wavelength (nm)')
ylabel('Spectrum (a. u.)')

figure(3);
subplot(211);
pcolor(t*1e15,z,abs(u1.^2)./max(max(abs(u1.^2))));
axis([t(1)*1e15 t(end)*1e15 0 z(end)]);
shading flat
xlabel('Time (fs)')
ylabel('Position (m)')
colorbar;
title('Temporal profile (a. u.)');
subplot(212);
pcolor(lambdav*1e9,z,10*log10(abs(uw1.^2)./max(max(abs(uw1.^2)))));
axis([100 1800 0 z(end)]);
clim([-30 0]);
colorbar;
shading flat
xlabel('Wavelength (nm)')
ylabel('Position (m)')
title('Spectrum (dB)');


[~,indextmax]=max(sum(sum(abs(uxtot.^2),1),2));
uxtot_attmax=squeeze(uxtot(:,:,indextmax));
uytot_attmax=squeeze(uytot(:,:,indextmax));

figure(4);
tiledlayout(2,2);
nexttile;
pcolor(t*1e15,x*1e6,squeeze(abs(uxtot(round(nx/2),:,:).^2)));
axis([t(1)*1e15 t(end)*1e15 -a*1e6 a*1e6]);
shading flat;
ylabel('Position (µm)');
xlabel('Time (fs)');
nexttile
plot(TETA2(R2<a),angle(uxtot_attmax(R2<a)),'o',TETA2(R2<a),angle(uytot_attmax(R2<a)),'+');
xlabel('Angle');
ylabel('Phase');
legend('x out','y out');
grid on;
nexttile
pcolor(x*1e6,x*1e6,squeeze(sum(abs(uxtot.^2),3)));
axis([-a*1e6 a*1e6 -a*1e6 a*1e6]);
axis square;
shading flat;
ylabel('Position (µm)');
xlabel('Position (µm)');
title('X polarization');
nexttile
pcolor(x*1e6,x*1e6,squeeze(sum(abs(uytot.^2),3)));
axis([-a*1e6 a*1e6 -a*1e6 a*1e6]);
axis square;
shading flat;
ylabel('Position (µm)');
xlabel('Position (µm)');
title('Y polarization');
    
% RDW plots
lambdaDW_start=100e-9;
lambdaDW_stop=400e-9;

uxtotDW=uxtot;
uxtotDWw=fftshift(fft(uxtotDW,[],3),3);
uxtotDWw(:,:,lambdav<lambdaDW_start | lambdav>lambdaDW_stop)=0;
uxtotDW=ifft(fftshift(uxtotDWw,3),[],3);

uytotDW=uxtot;
uytotDWw=fftshift(fft(uytotDW,[],3),3);
uytotDWw(:,:,lambdav<lambdaDW_start | lambdav>lambdaDW_stop)=0;
uytotDW=ifft(fftshift(uytotDWw,3),[],3);

EnergieDW=sum(sum(sum(abs(uxtotDW.^2)*dx*dx*dt)));

[~,indextmaxDW]=max(sum(sum(abs(uxtotDW.^2),1),2));
uxtot_attmaxDW=squeeze(uxtotDW(:,:,indextmaxDW));
uytot_attmaxDW=squeeze(uytotDW(:,:,indextmaxDW));


figure(5);
tiledlayout(3,2);
nexttile([1 2])
plot(t*1e15,squeeze(sum(sum(abs(uxtotDW.^2),1),2)*dx*dx)/max(squeeze(sum(sum(abs(uxtotDW.^2),1),2)*dx*dx)));
axis([-100 100 0 1.1]);
title('UV dispersive wave');
xlabel('Time (fs)');
ylabel('Intensity (a. u.)');
nexttile([1 2])
plot(lambdav(lambdav>0)*1e9,squeeze(sum(sum(abs(uxtotDWw(:,:,lambdav>0).^2),1),2)*dx*dx)/max(squeeze(sum(sum(abs(uxtotDWw(:,:,lambdav>0).^2),1),2)*dx*dx)));
axis([100 400 0 1.1])
xlabel('Wavelength (nm)')
ylabel('Spectrum (a. u.)')
nexttile
pcolor(x*1e6,x*1e6,squeeze(sum(abs(uxtotDW.^2),3)));
axis([-a*1e6 a*1e6 -a*1e6 a*1e6]);
axis square;
shading flat;
ylabel('Position (µm)');
xlabel('Position (µm)');
title('X polarization');
nexttile
plot(TETA2(R2<a),angle(uxtot_attmaxDW(R2<a)),'o',TETA2(R2<a),angle(uytot_attmaxDW(R2<a)),'+');
xlabel('Angle');
ylabel('Phase');
legend('x out','y out');
grid on;


% figure(10)
% pcolor(t*1e15,lambdav*1e9,spectrogram);shading flat;axis([-50 50 600 1500])
% xlabel('Time (fs)');
% ylabel('Wavelength (nm)');

LG102=R2/rOAM1.*exp(1i*TETA2).*exp(-R2.^2/rOAM1.^2);
LG202=(R2/rOAM1).^2.*exp(1i*2*TETA2).*exp(-R2.^2/rOAM2.^2);
InLG10DW=sum(sum(uxtot_attmaxDW.*conj(LG102)))./sqrt(sum(sum(abs(uxtot_attmaxDW.^2)))*sum(sum(abs(LG102.^2))));
InLG10=sum(sum(uxtot_attmax.*conj(LG102)))./sqrt(sum(sum(abs(uxtot_attmax.^2)))*sum(sum(abs(LG102.^2))));
InLG20DW=sum(sum(uxtot_attmaxDW.*conj(LG202)))./sqrt(sum(sum(abs(uxtot_attmaxDW.^2)))*sum(sum(abs(LG202.^2))));
InLG20=sum(sum(uxtot_attmax.*conj(LG202)))./sqrt(sum(sum(abs(uxtot_attmax.^2)))*sum(sum(abs(LG202.^2))));


disp(['Transmission = ' num2str(E_guided_v(end)/energie*100,'%.1f') ' %']);
disp(['Max B-integral per step Kerr= ' num2str(max(Bintv_Kerr),'%.3f') ' rad']);
disp(['Max B-integral per step Ion= ' num2str(max(Bintv_Ion),'%.3f') ' rad']);
disp(['Energy in DW = ' num2str(EnergieDW*1e6,'%.2f') ' µJ']);

switch input_type
    case 'LG10'
        disp(['Output at tmax in LG10 = ' num2str(abs(InLG10)*100,'%.1f') ' %']);
        disp(['Output at tmax in LG10 DW = ' num2str(abs(InLG10DW)*100,'%.1f') ' %']);
    case 'LG20'
        disp(['Output at tmax in LG20 = ' num2str(abs(InLG20)*100,'%.1f') ' %']);
        disp(['Output at tmax in LG20 DW = ' num2str(abs(InLG20DW)*100,'%.1f') ' %']);
end

% disp(['Ppeak/Pcrit_1bar = ' num2str(Pcrete/Pcrit_1bar,'%.2f') ]);
disp(['Compression GDD = ' num2str(phi2compopt*1e30,'%.0f') ' fs^2']);

toc;