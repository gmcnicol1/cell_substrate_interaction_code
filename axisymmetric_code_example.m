function axisymmetric_cell_code
clc
clear all
close all
format long
beep off

%% key parameters
% cell properties and geometry
R=12.5; % cell radius
A=pi*R^2; % cell adhered area
L=2*R; % length scale
H=2; % cell height
V=A*H;
N=50; % number of radial space steps
tmax=10000;
alpha=10; % ECM radius=alpha*cell radius

%dimensional position/time
Rn=5; % nucleus radius
r0=Rn;
r1=R;
rn=linspace(r0,r1,N+1)'; %dimensional positions
v0=1; % stall speed

Vcytoplasm=pi*(R^2-Rn^2)*H; %cytoplasm volume
Acytoplasm=pi*(R^2-Rn^2); %cytoplasm contact area with ECM

%dimensionless position/time
r0d=r0/L;
r1d=r1/L;
rd=linspace(r0d,r1d,N+1)'; % dimensionless position
hrd=rd(2)-rd(1); % dimensionless steps

%% dimensional biochemistry

% initial conditions
cG0=(100/Vcytoplasm)*ones(N+1,1);
cF0=zeros(N+1,1);
cSp0=zeros(N+1,1); 
cM0=(30/Vcytoplasm)*ones(N+1,1); 
cMp0=(0/Vcytoplasm)*ones(N+1,1); 
nf0=(100)*ones(N+1,1); 
nh0=zeros(N+1,1); 
nb0=zeros(N+1,1);
nA0=zeros(N+1,1);

% ligand patterning
NS0=1000;
ns0=NS0*ones(N+1,1);

% signalling proteins ICs
delta=4; % strength of increase in response to FA formation
cR0=(1e-3/Vcytoplasm)*ones(N+1,1); % ROCK
cRA0=(0/Vcytoplasm)*ones(N+1,1); % Activated ROCK 
cP0=(1/Vcytoplasm)*ones(N+1,1); % MLCP
cPP0=(0/Vcytoplasm)*ones(N+1,1); % MLCP-P
cC0=(1/Vcytoplasm)*ones(N+1,1); % cofilin
cCP0=(0/Vcytoplasm)*ones(N+1,1); % phosphorylated cofilin
cK0=(0.1/Vcytoplasm)*ones(N+1,1); % MLCK
cKP0=(0/Vcytoplasm)*ones(N+1,1); % MLCK-P 

% rate/diffusion constants

% structural/scaffolding proteins
Kpp=2e-3; %actin polymerisation
Kpm=1e-2; %actin depolymerisation
Kap=1e-2; %myosin activation
Kam=1e-2; %myosin inactivation
KFp=5e1; %adhesion maturation
KFm=1e-3; %focal adhesion disassembly

kmp=1e2; % myosin binding
kmm=1e-2; % myosin unbinding
khp=0.5; % affinity change
khm=5; % high affinity -> free integrin
kbp=1e-4; % nascent adhesion formation
kbm=1e-2; % nascent adhesion diassembly
DG=10; % actin monomer diffusivity
DF=0.03; % actin filament diffusivity
DM=1; % myosin II diffusivity
DMP=DM; % activated myosin II diffusivity
Dnf=0.1; % integrin diffusivity

% signaling proteins
KRp=1e-2; % ROCK activation
kRm=1e-1; % ROCK inactivation
K1p=1e-2;  % phosphorylation of MLCP
k1m=1e-2; % dephosphorylation of MLCP-P
K2p=1e-2; % phosphorylation of MLCK 
k2m=1e-2; % dephosphorylation of MLCK-P
K3p=5e-2; % phosphorylation of cofilin
k3m=1e-2; % dephosphorylation of phosphorylated cofilin
DP=15; % MLCP diffusivity
DPP=15; % MLCP-P diffusivity
DK=1; % MLCK diffusivity
DKP=1; % MLCK-P diffusivity
DC=10; % cofilin diffusivity
DCP=10; % phosphorylated cofilin diffusivity

%% scales
% concentration scaling
cA0bar=(2/(R^2-Rn^2))*trapz(rn,(cG0+cF0+cSp0).*rn); % actin
nI0bar=(2/(R^2-Rn^2))*trapz(rn,(nf0+nh0+nb0+nA0).*rn); % integrins/adhesions
nS0bar=(2/(R^2-Rn^2))*trapz(rn,ns0.*rn); % ligands
cMtot0bar=(2/(R^2-Rn^2))*trapz(rn,(cM0+cMp0).*rn); % myosin

cAtyp=cA0bar; % actin scale
nItyp=nI0bar; % integrin scale
cMtyp=cMtot0bar; % myosin scale
nStyp=nS0bar; % ligand scale
cPtyp=(2/(R^2-Rn^2))*trapz(rn,(cP0+cPP0).*rn); % myosin light chain phosphatase scale
cCtyp=(2/(R^2-Rn^2))*trapz(rn,(cC0+cCP0).*rn); % cofilin scale
cKtyp=(2/(R^2-Rn^2))*trapz(rn,(cK0+cKP0).*rn); % myosin light chain kinase scale
cRtyp=(2/(R^2-Rn^2))*trapz(rn,(cR0+cRA0).*rn); % ROCK scale

% time/velocity scaling
tchar=1/nStyp/kbp; % time scale
vchar=L/tchar; % velocity scale
v0d=v0*tchar/L; % dimensionless stall velocity

%% dimensionless biochemistry

% initial conditions
cG0d=cG0/cAtyp;
cF0d=cF0/cAtyp;
cSp0d=cSp0/cAtyp;
cM0d=cM0/cMtyp;
cMp0d=cMp0/cMtyp;
nf0d=nf0/nItyp;
nh0d=nh0/nItyp;
nb0d=nb0/nItyp;
nA0d=nA0/nItyp;
ns0d=ns0/nStyp;

cR0d=cR0/cRtyp;
cRA0d=cRA0/cRtyp;
cP0d=cP0/cPtyp;
cPP0d=cPP0/cPtyp;
cK0d=cK0/cKtyp;
cKP0d=cKP0/cKtyp;
cC0d=cC0/cCtyp;
cCP0d=cCP0/cCtyp;

% rate/diffusion constants
kmpMAt=kmp*cMtyp*cAtyp*tchar;
kmpAAt=kmp*cAtyp*cAtyp*tchar;
kmmt=kmm*tchar;
kmmAdMt=kmm*cAtyp*tchar/cMtyp;
khmt=khm*tchar;
kbpNt=kbp*tchar*nItyp;
kbpSt=kbp*tchar*nStyp;
kbmt=kbm*tchar;
kbmNdSt=kbm*nItyp*tchar/nStyp;
khpt=khp*tchar;
DGt=DG*tchar/L/L;
DFt=DF*tchar/L/L;
DMt=DM*tchar/L/L;
DMPt=DMP*tchar/L/L;
Dnft=Dnf*tchar/L/L;

kRmt=kRm*tchar;
k1mt=k1m*tchar;
k2mt=k2m*tchar;
k3mt=k3m*tchar;
DPt=DP*tchar/L/L;
DPPt=DPP*tchar/L/L;
DKt=DK*tchar/L/L;
DKPt=DKP*tchar/L/L;
DCt=DC*tchar/L/L;
DCPt=DCP*tchar/L/L;

%% Mechanics (dimensional and dimensionless)

% dimensional constants
ElF0=1000e-3; % filament stiffness
ElS0=5000e-3; % stress fibre stiffness 
ElP0=500e-3; % passive stiffness 
mu0=100e-3; % cytosol viscosity 
Elnuc=10000e-3; % nucleus stiffness

ElB0=0.01; % nascent adhesion stiffness constant
ElA0=0.05; % focal adhesion stiffness constant
bB0=1e-3; % nascent adhesion drag
bA0=5e-3; % focal adhesion drag

ElC0=1000; % ECM collagen stiffness
ElO0=0; % ECM 'other' stiffness
muE0=1e-3; % ECM viscosity

Fmax=1000e-3; % myosin II motor force

% dimensionless constants
ElS0d=ElS0/ElF0;
ElP0d=ElP0/ElF0;
Eldnuc=Elnuc/ElF0;
mu0d=mu0/ElF0/tchar;

ElB0d=ElB0*L*L/ElF0;
ElA0d=ElA0*L*L/ElF0;
bB0d=bB0*L*L/ElF0/tchar;
bA0d=bA0*L*L/ElF0/tchar;

ElC0d=ElC0/ElF0;
ElO0d=ElO0/ElF0;
muE0d=muE0/ElF0/tchar;

Fd=Fmax/ElF0;

% dimensionless functions;
fFd=cF0d; % filament stiffness dependence
fSd=cSp0d; % stress fibre stiffness dependence 
fPd=ones(N+1,1); % microtubules and nucleus
g1d=ones(N+1,1); % cytosol viscosity dependence
fCd=ones(N+1,1); % collagen
fOd=ones(N+1,1); % ECM 'other'
g2d=ones(N+1,1); % ECM viscosity dependence

fBd=nb0d; % nascent adhesions 
fAd=nA0d; % focal adhesions

hd=cSp0d; % active stress force

% full dimensionless functions
Eldcell=fFd+ElS0d*fSd+ElP0d*fPd;
EldECM=ElC0d*fCd+ElO0d*fOd;
mueffd=mu0d*g1d+Fd*hd/v0d;
muECMd=muE0d*g2d;
Eldadh=ElB0d*fBd+ElA0d*fAd; 
bBadhd=bB0d*fBd+bA0d*fAd+1e-5*ones(N+1,1);

% displacements
ur0(:,1)=zeros(N+3,1);
wr0(:,1)=zeros(N+3,1);

ur0d(:,1)=ur0/L;
wr0d(:,1)=wr0/L;

% stress initialisation
sigmacell=zeros(N+1,1);
sigmaECM=zeros(N+1,1);
lambda=zeros(N+1,1);
kfa=lambda;

%% Actin treadmilling (dimensional and dimensionless)

U0=1e-2;

U1=1+tanh(20*(rd(:)-r1d))-tanh(20*(rd(:)-r0d));

figure(22)
plot(rd,U1,'r','LineWidth',1.25)
xlabel('$x$','Interpreter','Latex','FontSize',16)
ylabel('$U$','Interpreter','Latex','FontSize',16)

U0sd=U0/vchar;
U1d=U0sd*U1;

%% ODE solver

tspan=[0,250,500,1000,1500,2000,2500,5000,tmax]/tchar; % timesteps

y0d=[cG0d;cF0d;cSp0d;cM0d;cMp0d;nf0d;nh0d;nb0d;nA0d;ns0d;cR0d;cRA0d;cP0d;cPP0d;cC0d;cCP0d;cK0d;cKP0d;ur0d;wr0d]; % dimensionless ICs

dcRdt=zeros(N+1,1);
dcRAdt=zeros(N+1,1);
dcPdt=zeros(N+1,1);
dcPPdt=zeros(N+1,1);
dcKdt=zeros(N+1,1);
dcKPdt=zeros(N+1,1);
dcCdt=zeros(N+1,1);
dcCPdt=zeros(N+1,1);
dcGdt=zeros(N+1,1);
dcFdt=zeros(N+1,1);
dcSpdt=zeros(N+1,1);
dcMdt=zeros(N+1,1);
dcMpdt=zeros(N+1,1);
dnfdt=zeros(N+1,1);
dnhdt=zeros(N+1,1);
dnbdt=zeros(N+1,1);
dnAdt=zeros(N+1,1);
dnsdt=zeros(N+1,1);
durddt=zeros(N+3,1);
dwrddt=zeros(N+3,1);

urd=ur0d;
wrd=wr0d;

Tsol=[];
Ysol=[];
tic
for j=1:length(tspan)-1
    tlocal=[tspan(j),tspan(j+1)];

    options=odeset('RelTol',1e-3,'AbsTol',1e-6);
    
    [tsol,ysol]=ode15s(@(t,y) axisymmetric_differential_equations_solver(t,y,N),tlocal,y0d,options);
    toc
    
    Tsol=[Tsol;tsol];
    Ysol=[Ysol;ysol];
    y0d=ysol(end,:);
    
    % cytoplasm (radial) stress
    sigmacell(1)=mueffd(1).*(-3*durddt(2)+4*durddt(3)-durddt(4))/2/hrd+Eldcell(1).*(-3*urd(2)+4*urd(3)-urd(4))/2/hrd+Fd*hd(1);
    sigmacell(2:N)=mueffd(2:N).*(durddt(4:N+2)-durddt(2:N))/2/hrd+Eldcell(2:N).*(urd(4:N+2)-urd(2:N))/2/hrd+Fd*hd(2:N);
    sigmacell(N+1)=mueffd(N+1).*(3*durddt(N+2)-4*durddt(N+1)+durddt(N))/2/hrd+Eldcell(N+1).*(3*urd(N+2)-4*urd(N+1)+urd(N))/2/hrd+Fd*hd(N+1);
    Scell(j+1,:)=sigmacell(:);
    
    % ECM (radial) stress
    sigmaECM(1)=muECMd(1).*(-3*dwrddt(2)+4*dwrddt(3)-dwrddt(4))/2/hrd+EldECM(1).*(-3*wrd(2)+4*wrd(3)-wrd(4))/2/hrd;
    sigmaECM(2:N)=muECMd(2:N).*(dwrddt(4:N+2)-dwrddt(2:N))/2/hrd+EldECM(2:N).*(wrd(4:N+2)-wrd(2:N))/2/hrd;
    sigmaECM(N+1)=muECMd(N+1).*(3*dwrddt(N+2)-4*dwrddt(N+1)+dwrddt(N))/2/hrd+EldECM(N+1).*(3*wrd(N+2)-4*wrd(N+1)+wrd(N))/2/hrd;
    SECM(j+1,:)=sigmaECM(:);
    
    straincell(j+1,:)=100*[(-3*urd(2)+4*urd(3)-urd(4))/2/hrd;(urd(4:N+2)-urd(2:N))/2/hrd;(3*urd(N+2)-4*urd(N+1)+urd(N))/2/hrd]; % strain cell
    strainECM(j+1,:)=100*[(-3*wrd(2)+4*wrd(3)-wrd(4))/2/hrd;(wrd(4:N+2)-wrd(2:N))/2/hrd;(3*wrd(N+2)-4*wrd(N+1)+wrd(N))/2/hrd]; % strain ECM
    
    outdurddt(j+1,:)=durddt;
    outdwrddt(j+1,:)=dwrddt;
end

outData1=[Tsol, Ysol(:,1:end)];
dlmwrite('axisymmetric_ode_output.dat',outData1,'delimiter','\t');

outData2=Scell;
dlmwrite('stress_cell.dat',outData2,'delimiter','\t');

outData3=SECM;
dlmwrite('stress_ECM.dat',outData3,'delimiter','\t');

outData4=straincell;
dlmwrite('strain_cell.dat',outData4,'delimiter','\t');

outData5=strainECM;
dlmwrite('strain_ECM.dat',outData5,'delimiter','\t');

outData6=outdurddt;
dlmwrite('outdurddt.dat',outData6,'delimiter','\t');

outData7=outdwrddt;
dlmwrite('outdwrddt.dat',outData7,'delimiter','\t');

axisymmetric_graphic_production_example(N,rn,rd,r0d,r1d,tchar,tspan,R,Rn,H,Acytoplasm,Vcytoplasm,cAtyp,nItyp,cRtyp,cPtyp,cKtyp,cCtyp)

%% functions
    function fval=axisymmetric_differential_equations_solver(t,y,N)

        cGd=y(1:N+1);
        cFd=y(N+2:2*N+2);
        cSpd=y(2*N+3:3*N+3);
        cMd=y(3*N+4:4*N+4);
        cMpd=y(4*N+5:5*N+5);
        nfd=y(5*N+6:6*N+6);
        nhd=y(6*N+7:7*N+7);
        nbd=y(7*N+8:8*N+8);
        nAd=y(8*N+9:9*N+9);
        nsd=y(9*N+10:10*N+10);
        cRd=y(10*N+11:11*N+11);
        cRAd=y(11*N+12:12*N+12);
        cPd=y(12*N+13:13*N+13);
        cPPd=y(13*N+14:14*N+14);
        cCd=y(14*N+15:15*N+15);
        cCPd=y(15*N+16:16*N+16);
        cKd=y(16*N+17:17*N+17);
        cKPd=y(17*N+18:18*N+18);
        urd=y(18*N+19:19*N+21);
        wrd=y(19*N+22:20*N+24);

        %% signalling proteins    
        kRp=KRp*(nbd(:)+delta*nAd(:));
        kRpt=kRp*tchar;

        dcRdt(:)=-kRpt(:).*cRd(:)+kRmt*cRAd(:);
        dcRAdt(:)=kRpt(:).*cRd(:)-kRmt*cRAd(:);

        k1p=K1p*cRAd(:);
        k2p=K2p*cRAd(:);
        k3p=K3p*cRAd(:);

        k1pt=k1p*tchar;
        k2pt=k2p*tchar;
        k3pt=k3p*tchar;
      
        dcPdt(1)=-k1pt(1).*cPd(1)+k1mt*cPPd(1)+2*DPt*(cPd(2)-cPd(1))/hrd^2;
        dcPPdt(1)=k1pt(1).*cPd(1)-k1mt*cPPd(1)+2*DPPt*(cPPd(2)-cPPd(1))/hrd^2;
        dcPdt(2:N)=-k1pt(2:N).*cPd(2:N)+k1mt*cPPd(2:N)+DPt*(cPd(3:N+1)-2*cPd(2:N)+cPd(1:N-1))/hrd^2+DPt*(cPd(3:N+1)-cPd(1:N-1))./rd(2:N)/2/hrd;
        dcPPdt(2:N)=k1pt(2:N).*cPd(2:N)-k1mt*cPPd(2:N)+DPPt*(cPPd(3:N+1)-2*cPPd(2:N)+cPPd(1:N-1))/hrd^2+DPPt*(cPPd(3:N+1)-cPPd(1:N-1))./rd(2:N)/2/hrd;
        dcPdt(N+1)=-k1pt(N+1).*cPd(N+1)+k1mt*cPPd(N+1)+2*DPt*(cPd(N)-cPd(N+1))/hrd^2;
        dcPPdt(N+1)=k1pt(N+1).*cPd(N+1)-k1mt*cPPd(N+1)+2*DPPt*(cPPd(N)-cPPd(N+1))/hrd^2;

        dcKdt(1)=-k2pt(1).*cKd(1)+k2mt*cKPd(1)+2*DKt*(cKd(2)-cKd(1))/hrd^2;
        dcKPdt(1)=k2pt(1).*cKd(1)-k2mt*cKPd(1)+2*DKPt*(cKPd(2)-cKPd(1))/hrd^2;
        dcKdt(2:N)=-k2pt(2:N).*cKd(2:N)+k2mt*cKPd(2:N)+DKt*(cKd(3:N+1)-2*cKd(2:N)+cKd(1:N-1))/hrd^2+DKt*(cKd(3:N+1)-cKd(1:N-1))./rd(2:N)/2/hrd;
        dcKPdt(2:N)=k2pt(2:N).*cKd(2:N)-k2mt*cKPd(2:N)+DKPt*(cKPd(3:N+1)-2*cKPd(2:N)+cKPd(1:N-1))/hrd^2+DKPt*(cKPd(3:N+1)-cKPd(1:N-1))./rd(2:N)/2/hrd;
        dcKdt(N+1)=-k2pt(N+1).*cKd(N+1)+k2mt*cKPd(N+1)+2*DKt*(cKd(N)-cKd(N+1))/hrd^2;
        dcKPdt(N+1)=k2pt(N+1).*cKd(N+1)-k2mt*cKPd(N+1)+2*DKPt*(cKPd(N)-cKPd(N+1))/hrd^2;

        dcCdt(1)=-k3pt(1).*cCd(1)+k3mt*cCPd(1)+2*DCt*(cCd(2)-cCd(1))/hrd^2;
        dcCPdt(1)=k3pt(1).*cCd(1)-k3mt*cCPd(1)+2*DCPt*(cCPd(2)-cCPd(1))/hrd^2;
        dcCdt(2:N)=-k3pt(2:N).*cCd(2:N)+k3mt*cCPd(2:N)+DCt*(cCd(3:N+1)-2*cCd(2:N)+cCd(1:N-1))/hrd^2+DCt*(cCd(3:N+1)-cCd(1:N-1))./rd(2:N)/2/hrd;
        dcCPdt(2:N)=k3pt(2:N).*cCd(2:N)-k3mt*cCPd(2:N)+DCPt*(cCPd(3:N+1)-2*cCPd(2:N)+cCPd(1:N-1))/hrd^2+DCPt*(cCPd(3:N+1)-cCPd(1:N-1))./rd(2:N)/2/hrd;
        dcCdt(N+1)=-k3pt(N+1).*cCd(N+1)+k3mt*cCPd(N+1)+2*DCt*(cCd(N)-cCd(N+1))/hrd^2;
        dcCPdt(N+1)=k3pt(N+1).*cCd(N+1)-k3mt*cCPd(N+1)+2*DCPt*(cCPd(N)-cCPd(N+1))/hrd^2;

        %% scaffolding proteins
        
        % actin polymerisation
        kpp=Kpp*cRAd;
        kpm=Kpm*cCd;
        
        % myosin II activation
        kap=Kap*cKPd;
        kam=Kam*cPd;
        
        % focal adhesion formation
        lambda=sqrt((urd(2:N+2)-wrd(2:N+2)).^2);
        kfa=0.5*lambda(:).^2;

        kFp=KFp*kfa;
        kFm=KFm*ones(N+1,1);

        %dimensionless reaction rates
        kppNt=kpp(:)*nI0bar*tchar;
        kpmt=kpm(:)*tchar;
        kapt=kap(:)*tchar;
        kamt=kam(:)*tchar;
        kFpt=kFp(:)*tchar;
        kFmt=kFm(:)*tchar;
        
        dcGdt(1)=-kppNt(1).*cGd(1).*(nbd(1)+nAd(1))+kpmt(1).*cFd(1)+kmmt*cSpd(1)+2*DGt*(cGd(2)-cGd(1))/hrd^2;
        dcGdt(2:N)=-kppNt(2:N).*cGd(2:N).*(nbd(2:N)+nAd(2:N))+kpmt(2:N).*cFd(2:N)+kmmt*cSpd(2:N)+DGt*(cGd(3:N+1)-2*cGd(2:N)+cGd(1:N-1))/hrd^2+DGt*(cGd(3:N+1)-cGd(1:N-1))/2/hrd./rd(2:N);
        dcGdt(N+1)=-kppNt(N+1).*cGd(N+1).*(nbd(N+1)+nAd(N+1))+kpmt(N+1).*cFd(N+1)+kmmt*cSpd(N+1)+2*DGt*(cGd(N)-cGd(N+1))/hrd^2;
        
        zeta1=cFd(:).*U1d(:).*rd(:);

        dcFdt(1)=kppNt(1).*cGd(1).*(nbd(1)+nAd(1))-kpmt(1).*cFd(1)-kmpMAt*(cFd(1)*cFd(1)*cMpd(1)+cFd(1)*cSpd(1)*cMpd(1))...
            +2*DFt*(cFd(2)-cFd(1)-hrd*U1d(1)*cFd(1)/DFt)/hrd/hrd+U1d(1)*cFd(1)/rd(1)-cFd(1)*U1d(1)/rd(1)-U1d(1)*U1d(1)*cFd(1)/DFt-cFd(1)*(-3*U1d(1)+4*U1d(2)-U1d(3))/2/hrd;
        dcFdt(2:N)=kppNt(2:N).*cGd(2:N).*(nbd(2:N)+nAd(2:N))-kpmt(2:N).*cFd(2:N)-kmpMAt*(cFd(2:N).*cFd(2:N).*cMpd(2:N)+cFd(2:N).*cSpd(2:N).*cMpd(2:N))+DFt*(cFd(3:N+1)-2*cFd(2:N)+cFd(1:N-1))/hrd^2+DFt*(cFd(3:N+1)-cFd(1:N-1))./rd(2:N)/2/hrd-(zeta1(3:N+1)-zeta1(1:N-1))/2/hrd./rd(2:N);
        dcFdt(N+1)=kppNt(N+1).*cGd(N+1).*(nbd(N+1)+nAd(N+1))-kpmt(N+1).*cFd(N+1)-kmpMAt*(cFd(N+1)*cFd(N+1)*cMpd(N+1)+cFd(N+1)*cSpd(N+1)*cMpd(N+1))+2*DFt*(cFd(N)-cFd(N+1)+hrd*U1d(N+1)*cFd(N+1)/DFt)/hrd^2+U1d(N+1)*cFd(N+1)/rd(N+1)-cFd(N+1)*U1d(N+1)/rd(N+1)-U1d(N+1)*U1d(N+1)*cFd(N+1)/DFt-cFd(N+1)*(3*U1d(N+1)-4*U1d(N)+U1d(N-1))/2/hrd;
           
        dcMdt(1)=-kapt(1)*cMd(1)+kamt(1)*cMpd(1)+2*DMt*(cMd(2)-cMd(1))/hrd^2;
        dcMdt(2:N)=-kapt(2:N).*cMd(2:N)+kamt(2:N).*cMpd(2:N)+DMt*(cMd(3:N+1)-2*cMd(2:N)+cMd(1:N-1))/hrd^2+DMt*(cMd(3:N+1)-cMd(1:N-1))./rd(2:N)/2/hrd;
        dcMdt(N+1)=-kapt(N+1)*cMd(N+1)+kamt(N+1).*cMpd(N+1)+2*DMt*(cMd(N)-cMd(N+1))/hrd^2;

        dcMpdt(1)=kapt(1)*cMd(1)-kamt(1)*cMpd(1)-kmpAAt*(cFd(1).*cFd(1).*cMpd(1)+cFd(1).*cSpd(1).*cMpd(1))+kmmAdMt*cSpd(1)+2*DMPt*(cMpd(2)-cMpd(1))/hrd^2;
        dcMpdt(2:N)=kapt(2:N).*cMd(2:N)-kamt(2:N).*cMpd(2:N)-kmpAAt*(cFd(2:N).*cFd(2:N).*cMpd(2:N)+cFd(2:N).*cSpd(2:N).*cMpd(2:N))+kmmAdMt*cSpd(2:N)+DMPt*(cMpd(3:N+1)-2*cMpd(2:N)+cMpd(1:N-1))/hrd^2+DMPt*(cMpd(3:N+1)-cMpd(1:N-1))./rd(2:N)/2/hrd;
        dcMpdt(N+1)=kapt(N+1)*cMd(N+1)-kamt(N+1).*cMpd(N+1)-kmpAAt*(cFd(N+1).*cFd(N+1).*cMpd(N+1)+cFd(N+1).*cSpd(N+1).*cMpd(N+1))+kmmAdMt*cSpd(N+1)+2*DMPt*(cMpd(N)-cMpd(N+1))/hrd^2;

        dnfdt(1)=-khpt*nfd(1)+khmt*nhd(1)+2*Dnft*(nfd(2)-nfd(1))/hrd^2;
        dnfdt(2:N)=-khpt*nfd(2:N)+khmt*nhd(2:N)+Dnft*(nfd(3:N+1)-2*nfd(2:N)+nfd(1:N-1))/hrd^2+Dnft*(nfd(3:N+1)-nfd(1:N-1))./rd(2:N)/2/hrd;
        dnfdt(N+1)=-khpt*nfd(N+1)+khmt*nhd(N+1)+2*Dnft*(nfd(N)-nfd(N+1))/hrd^2;

        dcSpdt(:)=kmpMAt*(cFd(:).*cFd(:).*cMpd(:)+cFd(:).*cSpd(:).*cMpd(:))-kmmt*cSpd(:);
        dnhdt(:)=khpt*nfd(:)-khmt*nhd(:)-kbpSt*nhd(:).*nsd(:)+kbmt*nbd(:);
        dnbdt(:)=kbpSt*nhd(:).*nsd(:)-kbmt*nbd(:)-kFpt(:).*nbd(:)+kFmt(:).*nAd(:);
        dnAdt(:)=kFpt(:).*nbd(:)-kFmt(:).*nAd(:);
        dnsdt(:)=-kbpNt*nhd(:).*nsd(:)+kbmNdSt*nbd(:);
        
        %% displacements

        fFd=cFd;
        fSd=cSpd;
        fBd=nbd;
        fAd=nAd;
        Eldcell=fFd+ElS0d*fSd+ElP0d*fPd;
        EldECM=ElC0d*fCd+ElO0d*fOd;
        Eldadh=ElB0d*fBd+ElA0d*fAd; 
        bBadhd=bB0d*fBd+bA0d*fAd+1e-5*ones(N+1,1);
        hd=cSpd;
        mueffd=mu0d*g1d+Fd*hd/v0d;
     
        gamma=0.05; % membrane tension constant    
        gammad=gamma/ElF0/L;
        
        V1r(1)=Eldnuc*urd(2)/rd(1)-Eldcell(1)*(urd(3)-urd(1))/2/hrd-Fd*hd(1);
        V1r(2)=(-(-3*Eldcell(1)+4*Eldcell(2)-Eldcell(3))/4/hrd/hrd-Eldcell(1)./rd(1)/2/hrd-Eldcell(1)/hrd/hrd).*urd(3)...
            +((-3*Eldcell(1)+4*Eldcell(2)-Eldcell(3))/4/hrd/hrd+Eldcell(1)./rd(1)/2/hrd-Eldcell(1)/hrd/hrd).*urd(1)...
            +(Eldcell(1)./rd(1)./rd(1)+2*Eldcell(1)/hrd/hrd+Eldadh(1)).*urd(2)...
            -Fd*(-3*hd(1)+4*hd(2)-hd(3))/2/hrd-Fd*hd(1)./rd(1)-Eldadh(1).*wrd(2);
        V1r(3:N+1)=(-(Eldcell(3:N+1)-Eldcell(1:N-1))/4/hrd/hrd-Eldcell(2:N)./rd(2:N)/2/hrd-Eldcell(2:N)/hrd/hrd).*urd(4:N+2)...
            +((Eldcell(3:N+1)-Eldcell(1:N-1))/4/hrd/hrd+Eldcell(2:N)./rd(2:N)/2/hrd-Eldcell(2:N)/hrd/hrd).*urd(2:N)...
            +(Eldcell(2:N)./rd(2:N)./rd(2:N)+2*Eldcell(2:N)/hrd/hrd+Eldadh(2:N)).*urd(3:N+1)...
            -Fd*(hd(3:N+1)-hd(1:N-1))/2/hrd-Fd*hd(2:N)./rd(2:N)-Eldadh(2:N).*wrd(3:N+1);
        V1r(N+2)=(-(3*Eldcell(N+1)-4*Eldcell(N)+Eldcell(N-1))/4/hrd/hrd-Eldcell(N+1)./rd(N+1)/2/hrd-Eldcell(N+1)/hrd/hrd).*urd(N+3)...
            +((3*Eldcell(N+1)-4*Eldcell(N)+Eldcell(N-1))/4/hrd/hrd+Eldcell(N+1)./rd(N+1)/2/hrd-Eldcell(N+1)/hrd/hrd).*urd(N+1)...
            +(Eldcell(N+1)./rd(N+1)./rd(N+1)+2*Eldcell(N+1)/hrd/hrd+Eldadh(N+1)).*urd(N+2)...
            -Fd*(3*hd(N+1)-4*hd(N)+hd(N-1))/2/hrd-Fd*hd(N+1)./rd(N+1)-Eldadh(N+1).*wrd(N+2);
        V1r(N+3)=-Eldcell(N+1)*(urd(N+3)-urd(N+1))/2/hrd-Fd*hd(N+1)+gammad*(1/(rd(N+1)+urd(N+1))-1/rd(N+1));

        V2r(1)=EldECM(1)*wrd(2)/rd(1)-EldECM(1)*(wrd(3)-wrd(1))/2/hrd;
        V2r(2)=(-(-3*EldECM(1)+4*EldECM(2)-EldECM(3))/4/hrd/hrd-EldECM(1)./rd(1)/2/hrd-EldECM(1)/hrd/hrd).*wrd(3)...
            +((-3*EldECM(1)+4*EldECM(2)-EldECM(3))/4/hrd/hrd+EldECM(1)./rd(1)/2/hrd-EldECM(1)/hrd/hrd).*wrd(1)...
            +(EldECM(1)./rd(1)./rd(1)+2*EldECM(1)/hrd/hrd+Eldadh(1)).*wrd(2)...
            -Eldadh(1).*urd(2);
        V2r(3:N+1)=(-(EldECM(3:N+1)-EldECM(1:N-1))/4/hrd/hrd-EldECM(2:N)./rd(2:N)/2/hrd-EldECM(2:N)/hrd/hrd).*wrd(4:N+2)...
            +((EldECM(3:N+1)-EldECM(1:N-1))/4/hrd/hrd+EldECM(2:N)./rd(2:N)/2/hrd-EldECM(2:N)/hrd/hrd).*wrd(2:N)...
            +(EldECM(2:N)./rd(2:N)./rd(2:N)+2*EldECM(2:N)/hrd/hrd+Eldadh(2:N)).*wrd(3:N+1)...
            -Eldadh(2:N).*urd(3:N+1);
        V2r(N+2)=(-(3*EldECM(N+1)-4*EldECM(N)+EldECM(N-1))/4/hrd/hrd-EldECM(N+1)./rd(N+1)/2/hrd-EldECM(N+1)/hrd/hrd).*wrd(N+3)...
            +((3*EldECM(N+1)-4*EldECM(N)+EldECM(N-1))/4/hrd/hrd+EldECM(N+1)./rd(N+1)/2/hrd-EldECM(N+1)/hrd/hrd).*wrd(N+1)...
            +(EldECM(N+1)./rd(N+1)./rd(N+1)+2*EldECM(N+1)/hrd/hrd+Eldadh(N+1)).*wrd(N+2)...
            -Eldadh(N+1).*urd(N+2);
        V2r(N+3)=(1+alpha*alpha)*EldECM(N+1)*wrd(N+2)/rd(N+1)/(1-alpha*alpha)-EldECM(N+1)*(wrd(N+3)-wrd(N+1))/2/hrd;
        
        % cell 
        M2r(1,1)=-mueffd(1)/2/hrd; 
        M2r(1,3)=mueffd(1)/2/hrd;
        M2r(2,1)=-(-3*mueffd(1)+4*mueffd(2)-mueffd(3))/4/hrd/hrd-mueffd(1)/rd(1)/2/hrd+mueffd(1)/hrd/hrd;
        M2r(2,2)=-mueffd(1)/rd(1)/rd(1)-2*mueffd(1)/hrd/hrd-bBadhd(1);
        M2r(2,3)=(-3*mueffd(1)+4*mueffd(2)-mueffd(3))/4/hrd/hrd+mueffd(1)/rd(1)/2/hrd+mueffd(1)/hrd/hrd;
        M2r(2,N+5)=bBadhd(1);
        % ECM
        M2r(N+4,N+4)=-muECMd(1)/2/hrd; % BC
        M2r(N+4,N+5)=-muECMd(1)/rd(1); % BC
        M2r(N+4,N+6)=muECMd(1)/2/hrd; % BC
        M2r(N+5,N+4)=-(-3*muECMd(1)+4*muECMd(2)-muECMd(3))/4/hrd/hrd-muECMd(1)/rd(1)/2/hrd+muECMd(1)/hrd/hrd;
        M2r(N+5,N+5)=-muECMd(1)/rd(1)/rd(1)-2*muECMd(1)/hrd/hrd-bBadhd(1);
        M2r(N+5,N+6)=(-3*muECMd(1)+4*muECMd(2)-muECMd(3))/4/hrd/hrd+muECMd(1)/rd(1)/2/hrd+muECMd(1)/hrd/hrd;
        M2r(N+5,2)=bBadhd(1);        
        for i=3:N+1
            % cell
            M2r(i,i-1)=-(mueffd(i)-mueffd(i-2))/4/hrd/hrd-mueffd(i-1)/rd(i-1)/2/hrd+mueffd(i-1)/hrd/hrd;
            M2r(i,i)=-mueffd(i-1)/rd(i-1)/rd(i-1)-2*mueffd(i-1)/hrd/hrd-bBadhd(i-1);
            M2r(i,i+1)=(mueffd(i)-mueffd(i-2))/4/hrd/hrd+mueffd(i-1)/rd(i-1)/2/hrd+mueffd(i-1)/hrd/hrd;
            M2r(i,i+N+3)=bBadhd(i-1);
            % ECM              
            M2r(N+3+i,N+2+i)=-(muECMd(i)-muECMd(i-2))/4/hrd/hrd-muECMd(i-1)/rd(i-1)/2/hrd+muECMd(i-1)/hrd/hrd;
            M2r(N+3+i,N+3+i)=-muECMd(i-1)/rd(i-1)/rd(i-1)-2*muECMd(i-1)/hrd/hrd-bBadhd(i-1);
            M2r(N+3+i,N+4+i)=(muECMd(i)-muECMd(i-2))/4/hrd/hrd+muECMd(i-1)/rd(i-1)/2/hrd+muECMd(i-1)/hrd/hrd;
            M2r(N+3+i,i)=bBadhd(i-1);
        end
        % cell 
        M2r(N+2,N+1)=-(3*mueffd(N+1)-4*mueffd(N)+mueffd(N-1))/4/hrd/hrd-mueffd(N+1)/rd(N+1)/2/hrd+mueffd(N+1)/hrd/hrd;
        M2r(N+2,N+2)=-mueffd(N+1)/rd(N+1)/rd(N+1)-2*mueffd(N+1)/hrd/hrd-bBadhd(N+1);
        M2r(N+2,N+3)=(3*mueffd(N+1)-4*mueffd(N)+mueffd(N-1))/4/hrd/hrd+mueffd(N+1)/rd(N+1)/2/hrd+mueffd(N+1)/hrd/hrd;
        M2r(N+2,2*N+5)=bBadhd(N+1);
        M2r(N+3,N+1)=-mueffd(N+1)/2/hrd;
        M2r(N+3,N+3)=mueffd(N+1)/2/hrd;
        % ECM
        M2r(2*N+5,2*N+4)=-(3*muECMd(N+1)-4*muECMd(N)+muECMd(N-1))/4/hrd/hrd-muECMd(N+1)/rd(N+1)/2/hrd+muECMd(N+1)/hrd/hrd;
        M2r(2*N+5,2*N+5)=-muECMd(N+1)/rd(N+1)/rd(N+1)-2*muECMd(N+1)/hrd/hrd-bBadhd(N+1);
        M2r(2*N+5,2*N+6)=(3*muECMd(N+1)-4*muECMd(N)+muECMd(N-1))/4/hrd/hrd+muECMd(N+1)/rd(N+1)/2/hrd+muECMd(N+1)/hrd/hrd;
        M2r(2*N+5,N+2)=bBadhd(N+1);
        M2r(2*N+6,2*N+4)=-muECMd(N+1)/2/hrd;
        M2r(2*N+6,2*N+5)=-(1+alpha*alpha)*muECMd(N+1)/rd(N+1)/(1-alpha*alpha);
        M2r(2*N+6,2*N+6)=muECMd(N+1)/2/hrd;
       
        Vr=[V1r';V2r'];
        
        Vrsol=(M2r)\(Vr); % (M*[dy/dt]=V)
        
        durddt=Vrsol(1:N+3);
        dwrddt=Vrsol(N+4:2*N+6);

        fval=[dcGdt;dcFdt;dcSpdt;dcMdt;dcMpdt;dnfdt;dnhdt;dnbdt;dnAdt;dnsdt;dcRdt;dcRAdt;dcPdt;dcPPdt;dcCdt;dcCPdt;dcKdt;dcKPdt;durddt;dwrddt];


    end
end