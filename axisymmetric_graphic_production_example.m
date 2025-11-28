function fval=axisymmetric_graphic_example(N,rn,rd,r0d,r1d,tchar,tspan,R,Rn,H,Acytoplasm,Vcytoplasm,cAtyp,nItyp,cRtyp,cPtyp,cKtyp,cCtyp)

%% graphic production
close all
format long
beep off

cmap=jet(256); % Get standard jet
pastel_factor=0.5; % 0 = white, 1 = original jet
cmap=(1-pastel_factor)+pastel_factor*cmap; % Blend toward white

A=load('axisymmetric_ode_output.dat'); % load data
stress_cell=load('stress_cell.dat'); % cell stress
stress_ECM=load('stress_ECM.dat'); % ECM stress
strain_cell=load('strain_cell.dat'); % cell strain
strain_ECM=load('strain_ECM.dat'); % ECM strain
outdurddt=load('outdurddt.dat'); % dur/dt
outdwrddt=load('outdwrddt.dat'); % dwr/dt

plotr=[1,0.85,0.7,0.55,0.4]; % selected spatial locations
for i=1:length(plotr) 
    y=plotr(i);
    [rdiff,p]=min(abs(rd-y));
    plotloc(i)=p;  
end

timeplot=zeros(length(tspan),1);

time=A(:,1); % simulation time
outcGd=A(:,2:N+2); % actin monomers
outcFd=A(:,N+3:2*N+3); % actin filaments
outcSpd=A(:,2*N+4:3*N+4); % stress fibres
outcMd=A(:,3*N+5:4*N+5); % inactive myosin II
outcMpd=A(:,4*N+6:5*N+6); % myosin II
outnfd=A(:,5*N+7:6*N+7); % free integrins
outnhd=A(:,6*N+8:7*N+8); % high affinity integrins
outnbd=A(:,7*N+9:8*N+9); % bound integrins 
outnAd=A(:,8*N+10:9*N+10); % focal adhesions
outnsd=A(:,9*N+11:10*N+11); % ECM ligands 
outcRd=A(:,10*N+12:11*N+12); % inactive ROCK 
outcRAd=A(:,11*N+13:12*N+13); % activated ROCK
outcPd=A(:,12*N+14:13*N+14); % MLCP
outcPPd=A(:,13*N+15:14*N+15); % MLCP-P
outcCd=A(:,14*N+16:15*N+16); % cofilin
outcCPd=A(:,15*N+17:16*N+17); % phosphorylated cofilin
outcKd=A(:,16*N+18:17*N+18); % MLCK
outcKPd=A(:,17*N+19:18*N+19); % MLCK-P
outurd=A(:,18*N+20:19*N+22); % ur
outwrd=A(:,19*N+23:20*N+25); % wr

% plotting temporal dynamics of proteins
tinset=2500; % zoomed in main figure

newcolors=[0 0 0;0.83 0.14 0.14;1.00 0.54 0.00;0.47 0.25 0.80;0.25 0.80 0.54];
timeplottingcolors=[0.5 0.5 0.5;0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.6350 0.0780 0.1840; 0 0 0];

for i=1:length(tspan)
    y=tspan(i);
    [tdiff,p]=min(abs(time-y));
    timeplot(i)=p;  
end

cGtot=[];
cFtot=[];
cSptot=[];
nftot=[];
nhtot=[];
nbtot=[];
nAtot=[];
cRAtot=[];
cPPtot=[];
cKPtot=[];
cCPtot=[];

for k=1:length(time)
    cGtot(k)=2*pi*H*trapz(rn,outcGd(k,:).*rn')/Vcytoplasm; % total fraction actin monomers
    cFtot(k)=2*pi*H*trapz(rn,outcFd(k,:).*rn')/Vcytoplasm; % total fraction actin filaments
    cSptot(k)=2*pi*H*trapz(rn,outcSpd(k,:).*rn')/Vcytoplasm; % total fraction actin in stress fibres
    nftot(k)=2*pi*trapz(rn,outnfd(k,:).*rn')/Acytoplasm; % total fraction free integrins
    nhtot(k)=2*pi*trapz(rn,outnhd(k,:).*rn')/Acytoplasm; % total fraction high-affinity integrins 
    nbtot(k)=2*pi*trapz(rn,outnbd(k,:).*rn')/Acytoplasm; % total fraction bound integrins
    nAtot(k)=2*pi*trapz(rn,outnAd(k,:).*rn')/Acytoplasm; % total fraction integrins in focal adhesions 
    cRAtot(k)=2*pi*H*trapz(rn,outcRAd(k,:).*rn')/Vcytoplasm; % total fraction of ROCK activated
    cPPtot(k)=2*pi*H*trapz(rn,outcPPd(k,:).*rn')/Vcytoplasm; % total fraction of MLCP phosphorylated
    cKPtot(k)=2*pi*H*trapz(rn,outcKPd(k,:).*rn')/Vcytoplasm; % total fraction of MLCK phosphorylated
    cCPtot(k)=2*pi*H*trapz(rn,outcCPd(k,:).*rn')/Vcytoplasm; % total fraction of cofilin phosphorylated
end

%% cell radius as a function of time
figure(1)
plot(time*tchar,r1d+outurd(:,N+2),'LineWidth',2.25)
colororder(timeplottingcolors)
set(gca,'fontsize',24)
xlabel('$t$ (s)','Interpreter','Latex','FontSize',36)
ylabel('$r_{c}/d_{c}$','Interpreter','Latex','FontSize',36)
xlim([0 tinset])
title('(a)','Interpreter','Latex','FontSize',36)
hold on
for i=2:length(tspan)
    plot(tspan(i)*tchar,r1d+outurd(timeplot(i),N+2),'o','Color',timeplottingcolors(i-1,:),'Linewidth',2.5)
    hold on
end
hold on
axes('Position',[.62 .625 .25 .25])
box on
plot(time*tchar,r1d+outurd(:,N+2),'LineWidth',2.25)
hold on
for i=length(tspan)-2:length(tspan)
    plot(tspan(i)*tchar,r1d+outurd(timeplot(i),N+2),'o','Color',timeplottingcolors(i-1,:),'Linewidth',2.5)
    hold on
end
set(gca,'Fontsize',16)
fig=gcf;
fig.Units='inches';
fig.Position=[0 0 8.5 6.5];        
fig.PaperUnits='inches';
fig.PaperSize=[8.5 6.5];          
fig.PaperPosition=[0 0 8.5 6.5];

%% global dynamics of cell biochemistry
figure(2)
pbaspect([3 1 1])
plot(time*tchar,100*nbtot,'Color',[0.3 0.25 0.93],'LineWidth',2.25)
hold on
plot(time*tchar,100*nAtot,'Color',[0.83 0.14 0.14],'LineWidth',2.25)
hold on
plot(time*tchar,100*cRAtot,'Color',[0.7 0 0.7],'LineWidth',2.25)
hold on
plot(time*tchar,100*cFtot,'Color',[1.00 0.54 0.00],'LineWidth',2.25)
hold on
plot(time*tchar,100*cSptot,'Color',[0.1 0.75 0.75],'LineWidth',2.25)
set(gca,'fontsize',24)
xlabel('$t$ (s)','Interpreter','Latex','FontSize',36)
ylabel('$f \: (\%)$','Interpreter','Latex','FontSize',36)
title('(b)','Interpreter','Latex','FontSize',36)
xlim([0 tinset])
ylim([0 100])
for i=2:length(tspan)
    plot(tspan(i)*tchar,100*nbtot(timeplot(i)),'o','Color',timeplottingcolors(i-1,:),'Linewidth',2.5)
    hold on
end
legend('Bound integrins','Focal adhesions','Activated ROCK','F-actin','VSFs','Interpreter','Latex','Location','northwest','Fontsize',16)
axes('Position',[.62 .615 .25 .25])
box on
plot(time*tchar,100*nbtot,'Color',[0.3 0.25 0.93],'LineWidth',2.25)
hold on
plot(time*tchar,100*nAtot,'Color',[0.83 0.14 0.14],'LineWidth',2.25)
hold on
plot(time*tchar,100*cRAtot,'Color',[0.7 0 0.7],'LineWidth',2.25)
hold on
plot(time*tchar,100*cFtot,'Color',[1.00 0.54 0.00],'LineWidth',2.25)
hold on
plot(time*tchar,100*cSptot,'Color',[0.1 0.75 0.75],'LineWidth',2.25)
colororder(newcolors)
hold on
for i=length(tspan)-2:length(tspan)
    plot(tspan(i)*tchar,100*nbtot(timeplot(i)),'o','Color',timeplottingcolors(i-1,:),'Linewidth',2.5)
    hold on
end
set(gca,'Fontsize',16)
fig=gcf;
fig.Units='inches';
fig.Position=[0 0 8.5 6.5];        
fig.PaperUnits='inches';
fig.PaperSize=[8.5 6.5];          
fig.PaperPosition=[0 0 8.5 6.5];


%% focal adhesion dynamics
timeplottingcolors=[0.5 0.5 0.5;0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.6350 0.0780 0.1840; 0 0 0];
figure(3)
plot(rd,outnAd(timeplot(2:end),:),'LineWidth',2.25)
set(gca,'fontsize',24)
colororder(timeplottingcolors)
hold on
xlabel('$r/d_{c}$','Interpreter','Latex','FontSize',36)
ylabel('$n_{A}/N_{I}$','Interpreter','Latex','FontSize',36)
title('(c)','Interpreter','Latex','FontSize',36)
ylim([0 inf])
hold on
hold on

axes('Position',[.2 .575 .25 .25])
box on
nangular=N;
theta=linspace(0,2*pi,nangular);
nAdend=outnAd(end,:);
Z=repmat(nAdend',1,nangular);
[R,Theta]=meshgrid(rd,theta);
X=R.*cos(Theta);
Y=R.*sin(Theta);
surf(X,Y,Z','EdgeColor','none')
surf(X,Y,Z','EdgeColor','none')
xlim([-0.55 0.55])
ylim([-0.55 0.55])
set(gca,'XTick',[-0.5 0 0.5])
set(gca,'YTick',[-0.5 0 0.5])
set(gca,'Fontsize',16)
view(2);
shading interp;
colormap(cmap);
colorbar;
box on
hold on
r_circle=rd(1);
x_circle=r_circle*cos(theta);
y_circle=r_circle*sin(theta);
fill(x_circle,y_circle,[0.35 0.35 0.35])
axis equal
hold on
plot(x_circle, y_circle, 'k--', 'LineWidth', 3)  
axis equal
fig=gcf;
fig.Units='inches';
fig.Position=[0 0 8.5 6.5];        
fig.PaperUnits='inches';
fig.PaperSize=[8.5 6.5];          
fig.PaperPosition=[0 0 8.5 6.5];

%% ROCK dynamics
figure(4)
plot(rd,outcRd(timeplot(2:end),:),'LineWidth',2.25)
set(gca,'fontsize',24)
colororder(timeplottingcolors)
hold on
plot(rd,outcRAd(timeplot(2:end),:),'--','LineWidth',2.25)
colororder(timeplottingcolors)
xlabel('$r/d_{c}$','Interpreter','Latex','FontSize',36)
ylabel('$c_{R}/C_{R},\:c_{R}^{+}/C_{R}$','Interpreter','Latex','FontSize',36)
title('(d)','Interpreter','Latex','FontSize',36)
hold on

axes('Position',[.2 .4 .25 .25])
box on
nangular=N;
theta=linspace(0,2*pi,nangular);
cRAdend=outcRAd(end,:);
Z=repmat(cRAdend',1,nangular);
[R,Theta]=meshgrid(rd,theta);
X=R.*cos(Theta);
Y=R.*sin(Theta);
surf(X,Y,Z','EdgeColor','none')
xlim([-0.55 0.55])
ylim([-0.55 0.55])
set(gca,'XTick',[-0.5 0 0.5])
set(gca,'YTick',[-0.5 0 0.5])
set(gca,'Fontsize',16)
view(2);
shading interp;
colormap(cmap);
colorbar;
box on
hold on
r_circle=rd(1);
x_circle=r_circle*cos(theta);
y_circle=r_circle*sin(theta);
fill(x_circle,y_circle,[0.35 0.35 0.35])
axis equal
hold on
plot(x_circle, y_circle, 'k--', 'LineWidth', 3)   % black dashed outline
axis equal
fig=gcf;
fig.Units='inches';
fig.Position=[0 0 8.5 6.5];        
fig.PaperUnits='inches';
fig.PaperSize=[8.5 6.5];          
fig.PaperPosition=[0 0 8.5 6.5];

%% actin filament dynamics
figure(5)
plot(rd,outcFd(timeplot(2:end),:),'LineWidth',2.25)
set(gca,'fontsize',24)
colororder(timeplottingcolors)
xlabel('$r/d_{c}$','Interpreter','Latex','FontSize',36)
ylabel('$c_{F}/C_{A}$','Interpreter','Latex','FontSize',36)
ylim([0 inf])
title('(e)','Interpreter','Latex','FontSize',36)
ylim([0 0.9])
hold on

axes('Position',[.2 .6 .25 .25])
box on
nangular=N;
theta=linspace(0,2*pi,nangular);
cFdend=outcFd(end,:);
Z=repmat(cFdend',1,nangular);
[R,Theta]=meshgrid(rd,theta);
X=R.*cos(Theta);
Y=R.*sin(Theta);
surf(X,Y,Z','EdgeColor','none')
surf(X,Y,Z','EdgeColor','none')
xlim([-0.55 0.55])
ylim([-0.55 0.55])
set(gca,'XTick',[-0.5 0 0.5])
set(gca,'YTick',[-0.5 0 0.5])
set(gca,'Fontsize',16)
view(2);
shading interp;
colormap(cmap);
colorbar;
box on
hold on
r_circle=rd(1);
x_circle=r_circle*cos(theta);
y_circle=r_circle*sin(theta);
fill(x_circle,y_circle,[0.35 0.35 0.35])
axis equal
hold on
plot(x_circle, y_circle, 'k--', 'LineWidth', 3)   % black dashed outline
axis equal
fig=gcf;
fig.Units='inches';
fig.Position=[0 0 8.5 6.5];        
fig.PaperUnits='inches';
fig.PaperSize=[8.5 6.5];          
fig.PaperPosition=[0 0 8.5 6.5];

%% stress fibre dynamics
figure(6)
plot(rd,outcSpd(timeplot(2:end),:),'LineWidth',2.25)
set(gca,'fontsize',24)
colororder(timeplottingcolors)
xlabel('$r/d_{c}$','Interpreter','Latex','FontSize',36)
ylabel('$c_{S}^{+}/C_{A}$','Interpreter','Latex','FontSize',36)
title('(f)','Interpreter','Latex','FontSize',36)
lg = legend('$t=250$ s','$t=500$ s','$t=1000$ s','$t=1500$ s',...
            '$t=2000$ s','$t=2500$ s','$t=5000$ s','$t=10000$ s',...
            'Interpreter','latex','NumColumns',3,'Location','northwest',...
            'FontSize',16);
lg.AutoUpdate = 'off';   
hold on

axes('Position',[.2 .325 .25 .25])
box on
nangular=N;
theta=linspace(0,2*pi,nangular);
cSpdend=outcSpd(end,:);
Z=repmat(cSpdend',1,nangular);
[R,Theta]=meshgrid(rd,theta);
X=R.*cos(Theta);
Y=R.*sin(Theta);
surf(X,Y,Z','EdgeColor','none')
surf(X,Y,Z','EdgeColor','none')
xlim([-0.55 0.55])
ylim([-0.55 0.55])
set(gca,'XTick',[-0.5 0 0.5])
set(gca,'YTick',[-0.5 0 0.5])
set(gca,'Fontsize',16)
view(2);
shading interp;
colormap(cmap);
colorbar;
box on
hold on
r_circle=rd(1);
x_circle=r_circle*cos(theta);
y_circle=r_circle*sin(theta);
fill(x_circle,y_circle,[0.35 0.35 0.35])
axis equal
hold on
plot(x_circle, y_circle, 'k--', 'LineWidth', 3)   % black dashed outline
axis equal
fig=gcf;
fig.Units='inches';
fig.Position=[0 0 8.5 6.5];        
fig.PaperUnits='inches';
fig.PaperSize=[8.5 6.5];          
fig.PaperPosition=[0 0 8.5 6.5];

%% cell displacement
figure(11)
plot(rd,outurd(timeplot(2:end),2:N+2),'LineWidth',2.25)
set(gca,'fontsize',24)
colororder(timeplottingcolors)
xlabel('$r/d_{c}$','Interpreter','Latex','FontSize',36)
ylabel('$u_{r}/d_{c}$','Interpreter','Latex','FontSize',36)
title('(a)','Interpreter','Latex','FontSize',36)
ylim([-inf 7.5e-3])
hold on

axes('Position',[.7 .67 .18 .18])
box on
plot(rd,strain_cell(2:end,:),'LineWidth',2.25)
set(gca,'Fontsize',16)
colororder(timeplottingcolors)
hold on
ylabel('$\partial{u_{r}}/\partial r \: (\%)$','Interpreter','Latex','FontSize',24)
axes('Position',[.2 .23 .25 .25])
box on
nangular=N;
theta=linspace(0,2*pi,nangular);
urdend=outurd(end,2:end-1);
Z=repmat(urdend',1,nangular);
[R,Theta]=meshgrid(rd,theta);
X=R.*cos(Theta);
Y=R.*sin(Theta);
surf(X,Y,Z','EdgeColor','none')
surf(X,Y,Z','EdgeColor','none')
xlim([-0.55 0.55])
ylim([-0.55 0.55])
set(gca,'XTick',[-0.5 0 0.5])
set(gca,'YTick',[-0.5 0 0.5])
set(gca,'Fontsize',16)
view(2);
shading interp;
colormap(cmap);
colorbar;
hold on
r_circle=rd(1);
x_circle=r_circle*cos(theta);
y_circle=r_circle*sin(theta);
Rn=rd(1);            
ur0c=urdend(1);          
nr_inner=80;            
r_inner=linspace(0,Rn,nr_inner);
theta_inner=linspace(0,2*pi,nangular);
[R_in,Theta_in]=meshgrid(r_inner,theta_inner);
U_inner=(R_in/Rn)*ur0c;
X_in=R_in.*cos(Theta_in);
Y_in=R_in.*sin(Theta_in);
surf(X_in,Y_in,U_inner,'EdgeColor','none')
box on
hold on
th=linspace(0,2*pi,400);
xN=rd(1)*cos(th);
yN=rd(1)*sin(th);
zN=ones(size(th))*max(U_inner(:))*1.001;
axis equal
hold on
plot3(xN,yN,zN,'k--','LineWidth',2)
axis equal
fig=gcf;
fig.Units='inches';
fig.Position=[0 0 8.5 6.5];        
fig.PaperUnits='inches';
fig.PaperSize=[8.5 6.5];          
fig.PaperPosition=[0 0 8.5 6.5];

%% cell stress
figure(12)
plot(rd,stress_cell(2:end,:),'LineWidth',2.25)
set(gca,'fontsize',24)
colororder(timeplottingcolors)
xlabel('$r/d_{c}$','Interpreter','Latex','FontSize',36)
ylabel('$\sigma_{c,r}/E_{0}^{F}$','Interpreter','Latex','FontSize',36)
title('(c)','Interpreter','Latex','FontSize',36)
ylim([0 inf])
hold on

axes('Position',[.6 .625 .25 .25])
box on
nangular=N;
theta=linspace(0,2*pi,nangular);
[R,Theta]=meshgrid(rd, theta);
X=R.*cos(Theta);
Y=R.*sin(Theta);
stress_cell_end=stress_cell(end,:);    
Z=repmat(stress_cell_end, nangular, 1);
surf(X, Y, Z, 'EdgeColor', 'none')
xlim([-0.55 0.55])
ylim([-0.55 0.55])
set(gca,'XTick',[-0.5 0 0.5])
set(gca,'YTick',[-0.5 0 0.5])
set(gca,'FontSize',16)
view(2);
shading interp;
colormap(cmap);
colorbar;
hold on
Rn=rd(1); 
ur0c=urdend(1); 
nr_inner=80; 
r_inner=linspace(0,Rn,nr_inner);
theta_inner=linspace(0,2*pi,nangular);
[R_in, Theta_in] = meshgrid(r_inner, theta_inner);
X_in=R_in.*cos(Theta_in);
Y_in=R_in.*sin(Theta_in);
stress_inner=10*(1/Rn)*ur0c*ones(size(R_in)); % update depending on assumed nucleus stiffness
surf(X_in,Y_in,stress_inner,'EdgeColor','none')
box on
hold on
hold on
th =linspace(0, 2*pi, 400);
xN=rd(1)*cos(th);
yN=rd(1)*sin(th);
zN=ones(size(th))*max(stress_inner(:))*1.001;
plot3(xN,yN,zN,'k--','LineWidth',2)
axis equal
fig=gcf;
fig.Units='inches';
fig.Position=[0 0 8.5 6.5];        
fig.PaperUnits='inches';
fig.PaperSize=[8.5 6.5];          
fig.PaperPosition=[0 0 8.5 6.5];


%% ECM displacement
figure(13)
plot(rd,outwrd(timeplot(2:end),2:N+2),'LineWidth',2.25)
set(gca,'fontsize',24)
colororder(timeplottingcolors)
xlabel('$r/d_{c}$','Interpreter','Latex','FontSize',36)
ylabel('$w_{r}/d_{c}$','Interpreter','Latex','FontSize',36)
title('(b)','Interpreter','Latex','FontSize',36)
hold on

axes('Position',[.225 .24 .18 .18])
box on
plot(rd,strain_ECM(2:end,:),'LineWidth',2.25)
set(gca,'Fontsize',16)
colororder(timeplottingcolors)
hold on
ylabel('$\partial w_{r}/\partial r \: (\%)$','Interpreter','Latex','FontSize',24)
fig=gcf;
fig.Units='inches';
fig.Position=[0 0 8.5 6.5];        
fig.PaperUnits='inches';
fig.PaperSize=[8.5 6.5];          
fig.PaperPosition=[0 0 8.5 6.5];

%% ECM stress
figure(14)
plot(rd,stress_ECM(2:end,:),'LineWidth',2.25)
set(gca,'fontsize',24)
colororder(timeplottingcolors)
xlabel('$r/d_{c}$','Interpreter','Latex','FontSize',36)
ylabel('$\sigma_{E,r}/E_{0}^{F}$','Interpreter','Latex','FontSize',36)
title('(d)','Interpreter','Latex','FontSize',36)
legend('$t=250$ s','$t=500$ s','$t=1000$ s','$t=1500$ s','$t=2000$ s','$t=2500$ s','$t=5000$ s','$t=10000$ s','Interpreter','Latex','NumColumns',3,'Location','northeast','Fontsize',16)
hold on
fig=gcf;
fig.Units='inches';
fig.Position=[0 0 8.5 6.5];        
fig.PaperUnits='inches';
fig.PaperSize=[8.5 6.5];          
fig.PaperPosition=[0 0 8.5 6.5];

end