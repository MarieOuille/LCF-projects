% 1D model of gas-filled MPC with adaptive beam size computed through
% ABCD matrix model of the Kerr lens - 1 gas 2 bulk materials - 4 bulk plates

clearvars;
tic;

% Initial condition

lambda0=5e-6;
tFWHM=200e-15;
phi2=0e-30;
phi3=0e-45;
energie=2e-6;

lambda1=0.84e-6;
tFWHM1=20e-15;
phi2_1=100000e-30;
phi3_1=0;
energie1=1e-9;

M2=1.2;

input_type='Gauss'; %Gauss or Exp or TwoPulses
File='Speck.dat';

% Materials
mat1='air';
pressure_1=1;
mat2='ZnSe';
mat3='CaF2';

% MPC geometry
f1= 0.1;%0.1; %focal length=R/2
n_roundtrip=10;

Lres= 0.39; %0.39; %distance entre les 2 miroirs 0<Lres<2R (stab cond)
z_mat3_start1=0.050; %these numbers must increase
z_mat3_stop1=0.050;
z_mat2_start1=0.180;
z_mat2_stop1=0.184;
z_mat2_start2=0.220 ;
z_mat2_stop2=0.224;
z_mat3_start2=0.210;
z_mat3_stop2=0.210;

MR=0.98; %mirror reflectivity
phi2mirror=-0e-30; %mirror GDD

% MPC mirror file
% zob=xlsread('mirror_spec.xlsx');
% lambdadata=zob(:,1)*1e-9;
% phi2data=zob(:,2)*1e-30;
% reflectivitydata=zob(:,3);
% phi2mirror=interp1(lambdadata,phi2data,lambdav);
% phi2mirror(isnan(phi2_PC94))=0;
% MR=interp1(lambdadata,reflectivitydata,lambdav);

% Temporal and spatial steps
dt=10e-15; %initially: 10e-15
T=4e-12; %initially: 4e-12
dz_MPC_mat1=5e-3;
dz_MPC_mat2=0.01e-3;
dz_MPC_mat3=0.01e-3;


%Compression parameters
phi2negmax=-60000e-30;
phi2posmax=60000e-30;
dphi2=100e-30;


%_____________________________________________

% Constants
c0=299792458;
hbar=1.05457173e-34;
eps0=8.85e-12;
w0=2*pi*c0/lambda0;

% Time-frequency variables
nt=2^(ceil(log2(T/dt)));
T=dt*nt;
t=-(nt/2-1)*dt:dt:(nt/2)*dt;
w=2*pi*(1/dt)*(-nt/2:nt/2-1)./nt;
dw=w(2)-w(1);
lambdav=c0*2*pi./(w+w0);
k0=2*pi/lambda0;
lmicron=lambdav*1e6;
lnano=lambdav*1e9;

% Material 1 (gas)

switch mat1
    case 'air'
        n2_1=pressure_1*3.6e-23;%air
        n_1=1+pressure_1*(0.05792105./(238.0185-lmicron.^(-2))+0.001679176./(57.362-lmicron.^(-2))); %air
        tau1_1=62.5e-15; %N2
        tau2_1=120e-15;  %N2
        fr_1=0.75;    %fr_1=0.75 for N2

    case 'N2'
        n2_1=pressure_1*4e-23; %N2 global
        n_1=1+pressure_1*(6.8552e-5+3.243157e-2./(144-lmicron.^-2)); %N2
        tau1_1=62.5e-15; %N2
        tau2_1=120e-15;  %N2
        fr_1=0.75;    %fr_1=0.75 for N2

    case 'Ar'
        n2_1=pressure_1*0.975e-23; %Ar
        n_1=1+pressure_1*(2.50141e-3./(91.012-lmicron.^-2)+5.00283e-4./(87.892-lmicron.^-2)+5.22343e-2./(214.02-lmicron.^-2)); %Ar
        tau1_1=1e-15;
        tau2_1=1e-15;
        fr_1=0;

    case 'Xe'
        n2_1=pressure_1*5.2e-23; %Xe
        n_1=1+pressure_1*(0.00322869./(46.301-lmicron.^-2)+0.00355393./(59.578-lmicron.^-2)+0.0606764./(112.74-lmicron.^-2)); %Xe
        tau1_1=1e-15;
        tau2_1=1e-15;
        fr_1=0;

    case 'Ne'
        n2_1=pressure_1*8.7e-25; %Ne
        n_1=1+pressure_1*(0.00128145./(184.661-lmicron.^-2)+0.0220486./(376.840-lmicron.^-2)); %Ne
        tau1_1=1e-15;
        tau2_1=1e-15;
        fr_1=0;

end

n_1=real(n_1);
n0_1=n_1(nt/2);
beta_1=(w+w0).*n_1/c0;
beta1_1=gradient(beta_1,w);
vg0_1=1/beta1_1(nt/2);

hr_1=(tau1_1^2+tau2_1^2)/(tau1_1*tau2_1^2).*exp(-t/tau2_1).*sin(t/tau1_1);
hr_1(t<0)=0;
hr_1=hr_1/sum(hr_1);
Resp_1=fr_1*hr_1;
Resp_1(t==0)=Resp_1(t==0)+(1-fr_1);

switch mat2
    case 'FS'
        n2_2=2.8e-20;            % FS
        n_2 = sqrt(1 + 0.6961663*lmicron.^2 ./ (lmicron.^2 - 0.0684043^2) + 0.4079426*lmicron.^2./(lmicron.^2 - 0.1162414^2) + 0.8974794*lmicron.^2./(lmicron.^2 - 9.896161^2)); %FS
        tau1_2=12.2e-15; %FS
        tau2_2=32e-15; %FS
        fr_2=0.18;    % FS

    case 'YAG'
        n2_2=8e-20;            % YAG
        n_2 = sqrt(1+2.28200./(1-0.01185./lmicron.^2)+3.27644./(1-282.734./lmicron.^2)); %YAG
        tau1_2=1e-15; %YAG
        tau2_2=1e-15; %YAG
        fr_2=0;    % YAG

    case 'Si'
        n2_2=4.5e-18;            % Si
        n_2 = sqrt(1+10.6684293./(1-(0.301516485./lmicron).^2)+0.0030434748./(1-(1.13475115./lmicron).^2)+1.54133408./(1-(1104./lmicron).^2));% %Si
        tau1_2=1/(2*pi*15.6e12); %Si
        tau2_2=1/(2*pi*105e9); %Si
        fr_2=0.026;    % Si

    case 'ZnSe'
        n2_2=12e-19;            % ZnSe
        n_2 = sqrt(1-0.689818+4.855169./(1-0.056359./lmicron.^2)+0.673922./(1-0.056336./lmicron.^2)+2.481890./(1-2222.114./lmicron.^2));%ZnSe
        tau1_2=1e-15; %ZnSe
        tau2_2=1e-15; %ZnSe
        fr_2=0;    % ZnSe

    case 'CaF2'
        n2_2=1.7e-20;            % CaF2
        n_2 = sqrt(1+0.33973+0.69913./(1-(0.09374./lmicron).^2)+0.11994./(1-(21.18./lmicron).^2)+4.35181./(1-(38.46./lmicron).^2)); %CaF2
        tau1_2=1e-15; %YAG
        tau2_2=1e-15; %YAG
        fr_2=0;    % YAG

    case 'Ge'
        n2_2=3e-17;            % Germanium
        n_2 = sqrt(1+0.4886331./(1-1.393959./lmicron.^2)+14.5142535./(1-0.1626427./lmicron.^2)+0.0091224./(1-752.190./lmicron.^2)); %Burnett et al 2016 (2-14um)
        tau1_2=1e-15; %YAG
        tau2_2=1e-15; %YAG
        fr_2=0;    % YAG

    case 'KBr'
        n2_2=4e-20;            % KBr
        n_2 = sqrt(1+0.39408+0.79221./(1-(0.146./lmicron).^2)+0.01981./(1-(0.173./lmicron).^2)+0.15587./(1-(0.187./lmicron).^2)+0.17673./(1-(60.61./lmicron).^2)+2.06217./(1-(87.72./lmicron).^2)); %BLi 1976 (0.2-42um)
        tau1_2=1e-15; %YAG
        tau2_2=1e-15; %YAG
        fr_2=0;    % YAG


end

n_2=real(n_2);
n0_2=n_2(nt/2);
beta_2=(w+w0).*n_2/c0;
beta1_2=gradient(beta_2,w);
vg0_2=1/beta1_2(nt/2);

hr_2=(tau1_2^2+tau2_2^2)/(tau1_2*tau2_2^2).*exp(-t/tau2_2).*sin(t/tau1_2);
hr_2(t<0)=0;
hr_2=hr_2/sum(hr_2);
Resp_2=fr_2*hr_2;
Resp_2(t==0)=Resp_2(t==0)+(1-fr_2);

switch mat3
    case 'CaF2'
        n2_3=1.7e-20;            % CaF2
        n_3 = sqrt(1+0.33973+0.69913./(1-(0.09374./lmicron).^2)+0.11994./(1-(21.18./lmicron).^2)+4.35181./(1-(38.46./lmicron).^2)); %CaF2
        tau1_3=1e-15; %YAG
        tau2_3=1e-15; %YAG
        fr_3=0;    % YAG

    case 'FS'
        n2_3=2.8e-20;            % FS
        n_3 = sqrt(1 + 0.6961663*lmicron.^2 ./ (lmicron.^2 - 0.0684043^2) + 0.4079426*lmicron.^2./(lmicron.^2 - 0.1162414^2) + 0.8974794*lmicron.^2./(lmicron.^2 - 9.896161^2)); %FS
        tau1_3=12.2e-15; %FS
        tau2_3=32e-15; %FS
        fr_3=0.18;    % FS

    case 'YAG'
        n2_3=8e-20;            % YAG
        n_3 = sqrt(1+2.28200./(1-0.01185./lmicron.^2)+3.27644./(1-282.734./lmicron.^2)); %YAG
        tau1_3=1e-15; %YAG
        tau2_3=1e-15; %YAG
        fr_3=0;    % YAG

    case 'Si'
        n2_3=4.5e-18;            % Si
        n_3 = sqrt(1+10.6684293./(1-(0.301516485./lmicron).^2)+0.0030434748./(1-(1.13475115./lmicron).^2)+1.54133408./(1-(1104./lmicron).^2));% %Si
        tau1_3=1/(2*pi*15.6e12); %Si
        tau2_3=1/(2*pi*105e9); %Si
        fr_3=0.026;    % Si

    case 'ZnSe'
        n2_3=12e-19;            % ZnSe
        n_3 = sqrt(1-0.689818+4.855169./(1-0.056359./lmicron.^2)+0.673922./(1-0.056336./lmicron.^2)+2.481890./(1-2222.114./lmicron.^2));%ZnSe
        tau1_3=1e-15; %ZnSe
        tau2_3=1e-15; %ZnSe
        fr_3=0;    % ZnSe

    case 'Ge'
        n2_3=3e-17;            % Germanium
        n_3 = sqrt(1+0.4886331./(1-1.393959./lmicron.^2)+14.5142535./(1-0.1626427./lmicron.^2)+0.0091224./(1-752.190./lmicron.^2)); %Burnett et al 2016 (2-14um)
        tau1_3=1e-15; %YAG
        tau2_3=1e-15; %YAG
        fr_3=0;    % YAG

    case 'KBr'
        n2_3=4e-20;            % KBr
        n_3 = sqrt(1+0.39408+0.79221./(1-(0.146./lmicron).^2)+0.01981./(1-(0.173./lmicron).^2)+0.15587./(1-(0.187./lmicron).^2)+0.17673./(1-(60.61./lmicron).^2)+2.06217./(1-(87.72./lmicron).^2)); %BLi 1976 (0.2-42um)
        tau1_3=1e-15; %YAG
        tau2_3=1e-15; %YAG
        fr_3=0;    % YAG
    
    case 'BaF2'
        n2_3=2e-20;            % BaF2
        n_3 = sqrt(1+0.643356./(1-(0.057789./lmicron).^2)+0.506762./(1-(0.10968./lmicron).^2)+3.8261./(1-(46.3864./lmicron).^2));
        tau1_3=1e-15; %YAG
        tau2_3=1e-15; %YAG
        fr_3=0;    % YAG
end

n_3=real(n_3);
n0_3=n_3(nt/2);
beta_3=(w+w0).*n_3/c0;
beta1_3=gradient(beta_3,w);
vg0_3=1/beta1_3(nt/2);

hr_3=(tau1_3^2+tau2_3^2)/(tau1_3*tau2_3^2).*exp(-t/tau2_3).*sin(t/tau1_3);
hr_3(t<0)=0;
hr_3=hr_3/sum(hr_3);
Resp_3=fr_3*hr_3;
Resp_3(t==0)=Resp_3(t==0)+(1-fr_3);


% MPC structure

dz_MPC_1=z_mat3_start1/round(z_mat3_start1/dz_MPC_mat1);
dz_MPC_1(isnan(dz_MPC_1))=dz_MPC_mat1;
dz_MPC_2=(z_mat3_stop1-z_mat3_start1)/round((z_mat3_stop1-z_mat3_start1)/dz_MPC_mat3);
dz_MPC_2(isnan(dz_MPC_2))=dz_MPC_mat3;
dz_MPC_3=(z_mat2_start1-z_mat3_stop1)/round((z_mat2_start1-z_mat3_stop1)/dz_MPC_mat1);
dz_MPC_3(isnan(dz_MPC_3))=dz_MPC_mat1;
dz_MPC_4=(z_mat2_stop1-z_mat2_start1)/round((z_mat2_stop1-z_mat2_start1)/dz_MPC_mat2);
dz_MPC_4(isnan(dz_MPC_4))=dz_MPC_mat2;
dz_MPC_5=(z_mat2_start2-z_mat2_stop1)/round((z_mat2_start2-z_mat2_stop1)/dz_MPC_mat1);
dz_MPC_5(isnan(dz_MPC_5))=dz_MPC_mat1;
dz_MPC_6=(z_mat2_stop2-z_mat2_start2)/round((z_mat2_stop2-z_mat2_start2)/dz_MPC_mat2);
dz_MPC_6(isnan(dz_MPC_6))=dz_MPC_mat2;
dz_MPC_7=(z_mat3_start2-z_mat2_stop2)/round((z_mat3_start2-z_mat2_stop2)/dz_MPC_mat1);
dz_MPC_7(isnan(dz_MPC_7))=dz_MPC_mat1;
dz_MPC_8=(z_mat3_stop2-z_mat3_start2)/round((z_mat3_stop2-z_mat3_start2)/dz_MPC_mat3);
dz_MPC_8(isnan(dz_MPC_8))=dz_MPC_mat3;
dz_MPC_9=(Lres-z_mat3_stop2)/round((Lres-z_mat3_stop2)/dz_MPC_mat1);
dz_MPC_9(isnan(dz_MPC_9))=dz_MPC_mat1;


z_MPC1=dz_MPC_1:dz_MPC_1:z_mat3_start1;
z_MPC2=z_mat3_start1+dz_MPC_2:dz_MPC_2:z_mat3_stop1;
z_MPC3=z_mat3_stop1+dz_MPC_3:dz_MPC_3:z_mat2_start1;
z_MPC4=z_mat2_start1+dz_MPC_4:dz_MPC_4:z_mat2_stop1;
z_MPC5=z_mat2_stop1+dz_MPC_5:dz_MPC_5:z_mat2_start2;
z_MPC6=z_mat2_start2+dz_MPC_6:dz_MPC_6:z_mat2_stop2;
z_MPC7=z_mat2_stop2+dz_MPC_7:dz_MPC_7:z_mat3_start2;
z_MPC8=z_mat3_start2+dz_MPC_8:dz_MPC_8:z_mat3_stop2;
z_MPC9=z_mat3_stop2+dz_MPC_9:dz_MPC_9:Lres;

z_MPC=[z_MPC1 z_MPC2 z_MPC3 z_MPC4 z_MPC5 z_MPC6 z_MPC7 z_MPC8 z_MPC9];

z_mat_MPC=ones(1,length(z_MPC));
z_mat_MPC(z_MPC>z_mat3_start1)=3;
z_mat_MPC(z_MPC>z_mat3_stop1)=1;
z_mat_MPC(z_MPC>z_mat2_start1)=2;
z_mat_MPC(z_MPC>z_mat2_stop1)=1;
z_mat_MPC(z_MPC>z_mat2_start2)=2;
z_mat_MPC(z_MPC>z_mat2_stop2)=1;
z_mat_MPC(z_MPC>z_mat3_start2)=3;
z_mat_MPC(z_MPC>z_mat3_stop2)=1;

M_rt=[1 0;-1/f1 1]*[1 z_mat3_start1/n0_1;0 1]*[1 (z_mat3_stop1-z_mat3_start1)/n0_3;0 1]*[1 (z_mat2_start1-z_mat3_stop1)/n0_1;0 1]*[1 (z_mat2_stop1-z_mat2_start1)/n0_2;0 1]*[1 (z_mat2_start2-z_mat2_stop1)/n0_1;0 1]*[1 (z_mat2_stop2-z_mat2_start2)/n0_2;0 1]*[1 (z_mat3_start2-z_mat2_stop2)/n0_1;0 1]*[1 (z_mat3_stop2-z_mat3_start2)/n0_3;0 1]*[1 (Lres-z_mat3_stop2)/n0_1;0 1]*[1 0;-1/f1 1]*[1 (Lres-z_mat3_stop2)/n0_1;0 1]*[1 (z_mat3_stop2-z_mat3_start2)/n0_3;0 1]*[1 (z_mat3_start2-z_mat2_stop2)/n0_1;0 1]*[1 (z_mat2_stop2-z_mat2_start2)/n0_2;0 1]*[1 (z_mat2_start2-z_mat2_stop1)/n0_1;0 1]*[1 (z_mat2_stop1-z_mat2_start1)/n0_2;0 1]*[1 (z_mat2_start1-z_mat3_stop1)/n0_1;0 1]*[1 (z_mat3_stop1-z_mat3_start1)/n0_3;0 1]*[1 z_mat3_start1/n0_1;0 1];
qini=1/((M_rt(2,2)-M_rt(1,1))/2/M_rt(1,2)+1i/abs(M_rt(1,2))*sqrt(1-((M_rt(1,1)+M_rt(2,2))/2)^2));
wwini=sqrt(lambda0/pi/imag(1/qini));

z_mat=zeros(2*n_roundtrip,length(z_MPC));
z_mat(1:2:end,:)=repmat(z_mat_MPC,n_roundtrip,1);
z_mat(2:2:end,:)=repmat(fliplr(z_mat_MPC),n_roundtrip,1);
z_mat=reshape(z_mat.',1,2*n_roundtrip*length(z_MPC));

% Initial condition processing

switch input_type
    case 'Gauss'
        tgauss=tFWHM/(sqrt(2*log(2)));
        e=exp(-(t/tgauss).^2);
        ew=fftshift(fft(e));
        ew=ew.*exp(-1i.*(w.^2/2*phi2+w.^3/6*phi3));
        e=ifft(ifftshift(ew));
        energieint=sum(abs(e.^2)*dt);
        ecrete=max(abs(e))*sqrt(energie/energieint);
        Pcrete=ecrete^2;
        e=e*sqrt(energie/energieint);


    case 'TwoPulses'
        tgauss=tFWHM/(sqrt(2*log(2)));
        tgauss1=tFWHM1/(sqrt(2*log(2)));
        deltaw1=2*pi*c0/lambda1-w0;

        e=exp(-(t/tgauss).^2);
        ew=fftshift(fft(e));
        ew=ew.*exp(-1i.*(w.^2/2*phi2+w.^3/6*phi3));
        e=ifft(ifftshift(ew));

        e1=exp(-(t/tgauss1).^2).*exp(1i*deltaw1*t);
        ew1=fftshift(fft(e1));
        ew1=ew1.*exp(-1i.*((w-deltaw1).^2/2*phi2_1+(w-deltaw1).^3/6*phi3_1));
        e1=ifft(ifftshift(ew1));

        energieint=sum(abs(e.^2)*dt);
        ecrete=max(abs(e))*sqrt(energie/energieint);
        Pcrete=ecrete^2;
        e=e*sqrt(energie/energieint);

        energieint1=sum(abs(e1.^2)*dt);
        ecrete1=max(abs(e1))*sqrt(energie1/energieint1);
        Pcrete1=ecrete1^2;
        e1=e1*sqrt(energie1/energieint1);

        e=e+e1;

    case 'Exp'
        aaa=load(File);
        ldata=aaa(:,1)*1e-9;
        aaa(:,2) = aaa(:,2)./max(aaa(:,2));
        ewdata=sqrt(aaa(:,2)).*exp(1i*aaa(:,3)); %might need to change phase sign
        ewdata=circshift(ewdata,0);
        ewinterp=interp1(ldata,ewdata,lambdav);
        ewinterp(isnan(ewinterp))=0;
        ew=ewinterp;
        e=fftshift(ifft(ifftshift(ew)));
        energieint=sum(abs(e.^2)*dt);
        ecrete=max(abs(e))*sqrt(energie/energieint);
        Pcrete=ecrete^2;
        e=e*sqrt(energie/energieint);
end


Broundtrip_anal=8*pi*n2_1*Pcrete/lambda0^2*atan(sqrt(Lres/(4*f1-Lres)));
Pcrit_1=3.77*pi*lambda0^2/8/pi^2/n2_1;
Pcrit_2=3.77*pi*lambda0^2/8/pi^2/n2_2;
%Zsf_2=0.367*zr./sqrt((sqrt(Pcrete/Pcrit_2)-0.852)^2-0.0219);

% Processing mirror dispersion

if isscalar(phi2mirror)
    phi2mirror=phi2mirror*ones(1,nt);
end

phimirror=cumsum(cumsum(phi2mirror))*dw*dw;
phi1mirror=(phimirror(end)-phimirror(1))/(nt*dw);
phimirror=phimirror-phi1mirror*w;


% Longitudinal spatial variable

dz_MPC=diff(z_MPC);
dz_MPC=[dz_MPC(1) dz_MPC];

dz=zeros(2*n_roundtrip,length(z_MPC));
dz(1:2:end,:)=repmat(dz_MPC,n_roundtrip,1);
dz(2:2:end,:)=repmat(fliplr(dz_MPC),n_roundtrip,1);
dz=reshape(dz.',1,2*n_roundtrip*length(z_MPC));

z=cumsum(dz);
nz=length(z);

ew=fftshift(fft(e))+sqrt(hbar*(w0+w)*nt/dt).*exp(1i*2*pi*rand(1,nt)); % bruit 1 photon par mode pour effets spontanés
u=zeros(nz,nt);
uw=zeros(nz,nt);
Bintv=zeros(1,nz);

Fmax=2*energie/wwini.^2/pi/M2;
Imax=2*Pcrete/wwini.^2/pi/M2;

q=zeros(1,nz);
ww=zeros(1,nz);
fKerr=zeros(1,nz);

for k=1:nz

    waitbar(k/nz);

    if k==1
        uw(k,:)=ew;
        q(k)=qini;
    else
        uw(k,:)=fftshift(fft(u(k-1,:)));
        q(k)=q(k-1);
    end

    if z_mat(k)==1
        n2=n2_1;
        n0=n0_1;
        beta=beta_1;
        vg0=vg0_1;
        Resp=Resp_1;
    elseif z_mat(k)==2
        n2=n2_2;
        n0=n0_2;
        beta=beta_2;
        vg0=vg0_2;
        Resp=Resp_2;
    elseif z_mat(k)==3
        n2=n2_3;
        n0=n0_3;
        beta=beta_3;
        vg0=vg0_3;
        Resp=Resp_3;
    end

    uw(k,:)=uw(k,:).*exp(1i*dz(k)/2.*(-beta+w/vg0));
    u(k,:)=ifft(fftshift(uw(k,:)));

    fKerr(k)=0.47*M2^2*pi*(lambda0./pi./imag(1/q(k)))^2/n2/max(abs(u(k,:).^2))/dz(k); %fKerr max at pulse peak (adjusted using Pcrit for gaussian beam)
    A=1-dz(k)/2/n0/fKerr(k);
    B=dz(k)/n0*(1-dz(k)/4/fKerr(k));
    C=-1./fKerr(k);
    D=1-dz(k)/2/n0/fKerr(k);
    q(k)=(A*q(k)+B)/(C*q(k)+D);
    ww(k)=sqrt(lambda0/pi/imag(1/q(k)));
    gamma_NL=2*pi*n2/lambda0/pi/ww(k)^2/M2;

    zobiwan=circshift(fftshift(ifft(fft(abs(u(k,:).^2)).*fft(Resp))),1,2); %Kerr + Raman
    u(k,:)=u(k,:).*exp(-dz(k)*1i*gamma_NL*(zobiwan));
    ft=gradient(1i*gamma_NL*zobiwan.*u(k,:),dt);
    u(k,:)=u(k,:)+1i/w0*dz(k)*ft;
    Bintv(k)=max(dz(k)*gamma_NL*(zobiwan));

    %    u(k,:)=u(k,:).*exp(-dz(k)*1i*gamma_NL*abs(u(k,:).^2)); % Instantaneous Kerr only
    %    ft=gradient(+1i*gamma_NL*abs(u(k,:).^2).*u(k,:),dt);
    %    u(k,:)=u(k,:)+1i/w0*dz(k)*ft;
    %    Bintv(k)=max(dz(k)*gamma_NL*abs(u(k,:).^2));


    uw(k,:)=fftshift(fft(u(k,:)));
    uw(k,:)=uw(k,:).*exp(1i*dz(k)/2.*(-beta+w/vg0));
    %uw(k,:)=uw(k,:).*exp(-(w/(dw*nt*0.45)).^40);
    u(k,:)=ifft(fftshift(uw(k,:)));
    u(k,:)=u(k,:).*exp(-(t/(dt*nt*0.45)).^40);

    if 2*max(abs(u(k,:).^2))/ww(k).^2/pi/M2>Imax
        Imax=2*max(abs(u(k,:).^2))/ww(k).^2/pi/M2;
    end

    if abs(z(k)/Lres-round(z(k)/Lres))<min(dz)/10
        A=1;
        B=0;
        C=-1/f1;
        D=1;
        q(k)=(A*q(k)+B)/(C*q(k)+D);
        uw(k,:)=uw(k,:).*sqrt(MR).*exp(-1i*phimirror);
        u(k,:)=ifft(fftshift(uw(k,:)));

        if 2*sum(abs(u(k,:).^2))*dt/ww(k).^2/pi/M2>Fmax
            Fmax=2*sum(abs(u(k,:).^2))*dt/ww(k).^2/pi/M2;
        end
    end

end

%% 

% Compression
phi2compv=phi2negmax:dphi2:phi2posmax;
%phi2compv=phi2negmax:dphi2:0; %on peut inverser ... 0:dphi2:phi2posmax si disp pos necessaire par ex
maxi=0;
Ncompopt=1;

for k=1:length(phi2compv)

    phi2comp=phi2compv(k);
    phi3comp=0;
    ecompw=fftshift(fft(u(end,:)));
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
ecompw=fftshift(fft(u(end,:)));
ecompw=ecompw.*exp(-1i.*(w.^2/2*phi2compopt+w.^3/6*phi3compopt));
ecomp=ifft(ifftshift(ecompw));
% [aa,bb]=max(abs(ecomp));
% ecomp=circshift(ecomp,nt/2-bb);


figure(1)
subplot(211)
plot(z,sqrt(M2)*ww*1e3);
xlabel('Propagation distance (m)');
ylabel('Beam radius (mm)');
title('Caustic');
subplot(212)
plot(z,cumsum(Bintv));
xlabel('Propagation distance (m)');
ylabel('B-integral (rad)');
title('B-integral');

figure(2);
subplot(211);
plot(lambdav(lambdav>0)*1e9,abs(uw(1,lambdav>0).^2)./max(abs(uw(1,lambdav>0).^2)),lambdav(lambdav>0)*1e9,abs(uw(end,lambdav>0).^2)./max(abs(uw(end,lambdav>0).^2)));
xlabel('Wavelength (nm)');
ylabel('Spectrum (a. u.)');
axis([3000 8000 0 1.1]);
legend('Input','Output');
subplot(212);
plot(t*1e15,abs(u(1,:).^2)./max(abs(u(1,:).^2)),t*1e15,abs(u(end,:).^2)./max(abs(u(end,:).^2)),t*1e15,abs(ecomp.^2)./max(abs(ecomp.^2)));
xlabel('Time (fs)');
ylabel('Intensity (a. u.)');
axis([-1000 1000 0 1.1]);
legend('Input','Output','Compressed');

% figure(3);
% subplot(311)
% plot(z,sqrt(M2)*ww*1e3);
% xlabel('Propagation distance (m)');
% ylabel('Beam radius (mm)');
% %title('Caustic');
% subplot(312)
% plot(lambdav*1e9,abs(uw(end,:).^2)./max(abs(uw(end,:).^2)));
% xlabel('Wavelength (nm)');
% ylabel('Spectrum (a. u.)');
% axis([980 1080 0 1.1]);
% subplot(313);
% plot(t*1e12,abs(u(end,:).^2)./max(abs(u(end,:).^2)));
% xlabel('Time (ps)');
% ylabel('Intensity (a. u.)');
% axis([-2 2 0 1.1]);

% figure(4)
% plot(lambdav*1e9,10*log10(abs(uw(end,:).^2)./max(abs(uw(end,:).^2))));
% xlabel('Wavelength (nm)');
% ylabel('Spectrum (a. u.)');
% axis([800 1400 -80 0]);


toc;


disp(['Transmission = ' num2str(100*sum(abs(u(end,:).^2))./sum(abs(e.^2)),'%.1f') ' %']);
disp(['Ppeak/Pcrit1 = ' num2str(Pcrete/Pcrit_1,'%.2f') ]);
disp(['Ppeak/Pcrit2 = ' num2str(Pcrete/Pcrit_2,'%.2f') ]);
disp(['Compression GDD = ' num2str(phi2compopt*1e30,'%.0f') ' fs^2']);
disp(['Max B-integral per step = ' num2str(max(Bintv),'%.3f') ' rad']);
disp(['Max Fluence on mirror = ' num2str(Fmax*1e-4,'%.3f') ' J/cm^2']);
disp(['Max intensity = ' num2str(Imax*1e-4,'%.1e') ' W/cm^2']);



%
% %RAMAN
% % Isolating signals
%
% w0pump=0;
% BWpump=20/tFWHM;
% w1pump=w0pump-BWpump;
% w2pump=w0pump+BWpump;
% [~,indexpump1]=min(abs(w-w1pump));
% [~,indexpump2]=min(abs(w-w2pump));
% Epump=sum(abs(uw(:,indexpump1:indexpump2)).^2,2)*dt/nt;
% epumpw=uw;
% epumpw(:,1:indexpump1)=0;
% epumpw(:,indexpump2:end)=0;
% epump=ifft(fftshift(epumpw,2),[],2);
%
% w0Raman=-1/tau1_2;
% BWRaman=max(5/tau2_2,BWpump);
% w1Raman=+w0Raman-BWRaman;
% w2Raman=+w0Raman+BWRaman;
% [~,indexRaman1]=min(abs(w-w1Raman));
% [~,indexRaman2]=min(abs(w-w2Raman));
% ERaman=sum(abs(uw(:,min(indexRaman1,indexRaman2):max(indexRaman1,indexRaman2))).^2,2)*dt/nt;
% eramanw=uw;
% eramanw(:,1:indexRaman1)=0;
% eramanw(:,indexRaman2:end)=0;
% eraman=ifft(fftshift(eramanw,2),[],2);
%
% w0Raman2nd=-2/tau1_2;
% w1Raman2nd=+w0Raman2nd-BWRaman;
% w2Raman2nd=+w0Raman2nd+BWRaman;
% [~,indexRaman12nd]=min(abs(w-w1Raman2nd));
% [~,indexRaman22nd]=min(abs(w-w2Raman2nd));
% ERaman2nd=sum(abs(uw(:,indexRaman12nd:indexRaman22nd)).^2,2)*dt/nt;
% eramanw2nd=uw;
% eramanw2nd(:,1:indexRaman12nd)=0;
% eramanw2nd(:,indexRaman22nd:end)=0;
% eraman2nd=ifft(fftshift(eramanw2nd,2),[],2);
%
% w0Raman3rd=-3/tau1_2;
% w1Raman3rd=+w0Raman3rd-BWRaman;
% w2Raman3rd=+w0Raman3rd+BWRaman;
% [~,indexRaman13rd]=min(abs(w-w1Raman3rd));
% [~,indexRaman23rd]=min(abs(w-w2Raman3rd));
% ERaman3rd=sum(abs(uw(:,indexRaman13rd:indexRaman23rd)).^2,2)*dt/nt;
% eramanw3rd=uw;
% eramanw3rd(:,1:indexRaman13rd)=0;
% eramanw3rd(:,indexRaman23rd:end)=0;
% eraman3rd=ifft(fftshift(eramanw3rd,2),[],2);
%
% w0RamanAS=1/tau1_2;
% w1RamanAS=+w0RamanAS-BWRaman;
% w2RamanAS=+w0RamanAS+BWRaman;
% [~,indexRaman1AS]=min(abs(w-w1RamanAS));
% [~,indexRaman2AS]=min(abs(w-w2RamanAS));
% ERamanAS=sum(abs(uw(:,indexRaman1AS:indexRaman2AS)).^2,2)*dt/nt;
% eramanwAS=uw;
% eramanwAS(:,1:indexRaman1AS)=0;
% eramanwAS(:,indexRaman2AS:end)=0;
% eramanAS=ifft(fftshift(eramanwAS,2),[],2);
%
% %Compression 1st Stokes
% phi2compv=-phi2-10*phi2/100:phi2/100:-phi2+10*phi2/100;
% maxi=0;
% Ncompopt=1;
%
% for k=1:length(phi2compv)
%
% phi2comp=phi2compv(k);
% phi3comp=0;
% ecompw=fftshift(fft(eraman(end,:)));
% ecompw=ecompw.*exp(-1i.*((w-w0Raman).^2/2*phi2comp+(w-w0Raman).^3/6*phi3comp));
% ecomp=ifft(fftshift(ecompw));
% maximoum=max(abs(ecomp.^2));
%
% if maximoum>maxi
% Ncompopt=k;
% maxi=maximoum;
% end
%
% end
%
% phi2compopt=phi2compv(Ncompopt);
% phi3compopt=0;
% ecompw=fftshift(fft(eraman(end,:)));
% ecompw=ecompw.*exp(-1i.*(w.^2/2*phi2compopt+w.^3/6*phi3compopt));
% ecomp=ifft(fftshift(ecompw));
%
% %Plots Raman
%
% figure(4);
% subplot(311);
% plot(z,sum(abs(u.^2),2)/sum(abs(u(1,:).^2)));
% subplot(312);
% plot(z,ERaman,z,ERaman2nd,z,ERaman3rd,z,ERamanAS);%,z,ERaman(1).*exp(g0*Icrete*z));
% axis([0 z(end) 0 energie]);
% legend('1st Stokes','2nd Stokes','3rd Stokes','anti Stokes');
% xlabel('Propagation distance (m)');
% ylabel('Energy (J)');
% subplot(313);
% plot(z,Epump);
% axis([0 z(end) 0 energie]);
% xlabel('Propagation distance (m)');
% ylabel('Energy (J)');
%
% figure(5);
% subplot(211)
% plot(t*1e12,abs(u(1,:).^2),t*1e12,abs(epump(end,:).^2),t*1e12,abs(eraman(end,:).^2),t*1e12,abs(eraman2nd(end,:).^2),t*1e12,abs(eraman3rd(end,:).^2),t*1e12,abs(eramanAS(end,:).^2));
% legend('input','pump','1st Stokes','2nd Stokes','3rd Stokes','anti Stokes');
% xlabel('Time (ps)');
% ylabel('Power (W)');
% subplot(212)
% plot(lambdav*1e9,10*log10(abs(uw(1,:).^2)./max(abs(uw(1,:).^2))),lambdav*1e9,10*log10(abs(uw(end,:).^2)./max(abs(uw(1,:).^2))));
% axis([600 2500 -30 0]);
% xlabel('Wavelength (nm)');
% ylabel('Power (dB)');
%
% figure(6);
% subplot(211)
% plot(t*1e12,abs(eraman(end,:).^2),t*1e12,abs(ecomp(end,:).^2));
% legend('1st Stokes','1st Stokes compressed');
% xlabel('Time (ps)');
% ylabel('Power (W)');
% subplot(212)
% plot(lambdav*1e9,abs(eramanw(end,:).^2)./max(abs(eramanw(end,:).^2)));
% axis([lambdav(indexRaman2)*1e9 lambdav(indexRaman1)*1e9 0 1]);
% xlabel('Wavelength (nm)');
% ylabel('Power (u. a.)');


