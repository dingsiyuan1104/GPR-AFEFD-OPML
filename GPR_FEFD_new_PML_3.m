%     2D Finite Element Frequency Domain Simulation for GPR with PML 
%                the optimal PML & the general PML
% -------a single shot point----------a single frequency of 100MHz------- %  
% 
% 2D transverse magnetic(TM) mode waves
% d(dEz/dx)/dx+d(dEz/dy)/dy+k0^2*(epr-1i*sig/omega/eps0)*Ez=1i*omega*miu0*J;
%       k0^2 = omega^2*eps0*miu0;
%      d(A*du/dx)/dx+d(A*du/dy)/dy+C*k^2*u=0;
%      k^2 = omega^2*miu*eps0(epr-1i*sig/omega/eps0);
%      A = Sy/Sx ; B = Sx/Sy ; C = Sx*Sy
% genera PML
%      Sx = 1 - 1i*sig_x/omega/eps0;
%      sig_x = sigma_max*(x/d)^m;sigma_max=-(m+1)*ln(R)/(2*d*eta0*sqrt(eps));
% optimal PML
%      Sx = 1 - 1i*sig_x/k;
%      sig_x = 1/(d-x);
%
% -----------2018.10£¨Deshan Feng, Siyuan Ding, and Xun Wang£©----------- %
close all;
clc;
clear;
addpath Solvers
tic
%% Load model parameters (area size and boundary partition)
X = -4.2:0.05:4.2;
Y = -4.2:0.05:4.2;
ds = 0.05;          % Discrete mesh spacing
xl = -4; xr = 4;    % The left and right boundaries of the calculation region in the X direction
yl = -4; yr = 4;    % The left and right boundaries of the calculation region in the Y direction
fprintf('thickness of pml is:');
d =   0.2;            % thickness of pml
%% Generate the array of nodes
for i=1:length(Y);
    for j=1:length(X);
        n=(i-1)*length(X)+j;
        node(n,1)=X(j);
        node(n,2)=Y(i);
    end
end
%% Mesh generation
elem=delaunay(node);  
% coordinates of elements
xq=(node(elem(:,1),1)+node(elem(:,2),1)+node(elem(:,3),1))/3;
yq=(node(elem(:,1),2)+node(elem(:,2),2)+node(elem(:,3),2))/3;
xv = [xl,xr,xr,xl,xl];
yv = [yl,yl,yr,yr,yl];
[in1,~] = inpolygon(xq,yq,xv,yv);      % simulation region
%% Parameters of the medium
elem(:,4)=3;          % relative permittivity
elem(:,5)=0.00005;    % conductivity
%% Parameters of the antenna
freq=100.0*10^6;            % main frequency
ep0=1/(36*pi)*1e-9;
miu=4*pi*1e-7;
omega = 2 * pi *freq;  
ep=elem(:,4);
sig=elem(:,5);
k2 = omega^2.*ep0*miu*(ep-1i*sig/omega/ep0);
%% Location of shot point (a single shot point)
sitexy(:,1)=0;              % x_coordinate
sitexy(:,2)=0;              % y_coordinate
%% mesh parameters for AFEM
mesh = getmesh_AFEM_tri(node,elem,[],[]);              % node;elem;type;solu
%% Pre loading of excitation source
nsorc = size(sitexy,1);                       % %Number of the shot points
[~,site]=ismember(sitexy,mesh.node,'rows');   % Sorting of the node of the shot points 
ND = size(mesh.node,1);                       % Number of nodes
RHS = sparse(ND,nsorc);                       % Source vector
for i = 1:nsorc
    RHS(site(i),i) = 1i*omega*miu;           
end
%% Forward modeling with the general PML
m=4;                % degree of polynomial grading (usually between 3 and 4)
R=exp(-16);         % theoretical reflection coefficient of the target with normal incidence
pml = pml_freq_fem2d_tri(mesh.elem,mesh.node,m,R,xl,xr,yl,yr,d);     % Calculate sigx and sigy
forward=pre_freq_AFEM2d_tri(mesh,pml,freq);                          % Generate the structure 
K = forward.K;
u=zeros(ND,nsorc);
for i = 1:nsorc
    u(:,i)=opti_zmumps(K,RHS(:,i));           % solutions with different shot points
end
toc
%% Forward modeling with the optimal PML  
pml_new = pml_new_freq_fem2d_tri(mesh.elem,mesh.node,xl,xr,yl,yr,d); % Calculate sigx and sigy
forward_new=pre_freq_AFEM2d_tri_new_pml(mesh,pml_new,freq);          % Generate the structure 
K_new=forward_new.K;
u_new=zeros(ND,nsorc);
for i = 1:nsorc
     u_new(:,i)=opti_zmumps(K_new,RHS(:,i));  % solutions with different shot points
end
toc
%% analytical solution  
current = 1.0;
k = sqrt(mean(k2(:)));
rr = ((node(:,1)-node(site,1)).^2+(node(:,2)-node(site,2)).^2).^(1/2); 
Utrue(:,1) = -current .* omega * miu .*besselh(0,2, k * rr )/4;
%% ----------------------------------Wave field diagram-----------------------------------%%
figure(3)    
set(gcf,'Position',[12 53 854 869]);      %position of picture
%----------------------------------------Analytical solution----------------------------------------------%
subplot(3,3,1)   %------Real part------%
Utrue_xy=reshape(full(real(Utrue)),169,169);
imagesc(X,Y,Utrue_xy);
title('(a) Re(Ut)','FontSize',11,'fontweight','bold');
xlabel('x(m)');
ylabel('y(m)');
view(2);axis image;colormap jet;
caxis([-100,100]);
set(gca,'FontSize',11,'fontweight','bold');
subplot(3,3,2)   %------Imaginary part------%
Utrue_xy=reshape(full(imag(Utrue)),169,169);
imagesc(X,Y,Utrue_xy);
title('(b) Im(Ut)','FontSize',11,'fontweight','bold');
xlabel('x(m)');
ylabel('y(m)');
view(2);axis image;colormap jet;
caxis([-100,100]);
set(gca,'FontSize',11,'fontweight','bold');
subplot(3,3,3)   %------Absolute value------%
Utrue_xy=reshape(full(abs(Utrue)),169,169);
imagesc(X,Y,Utrue_xy);
title('(c) Abs(Ut)','FontSize',11,'fontweight','bold');
xlabel('x(m)');
ylabel('y(m)');
view(2);axis image;colormap jet;
caxis([-100,100]);
set(gca,'FontSize',11,'fontweight','bold');
%----------------------------------------General PML----------------------------------------------%
subplot(3,3,4)   %------Real part------%
u_xy=reshape(full(real(u(:))),169,169);
imagesc(X,Y,u_xy);
title('(d) Re(Ug)','FontSize',11,'fontweight','bold');
xlabel('x(m)');
ylabel('y(m)');
view(2);axis image;colormap jet;
caxis([-100,100]);
set(gca,'FontSize',11,'fontweight','bold');
subplot(3,3,5)   %------Imaginary part------%
u_xy=reshape(full(imag(u(:))),169,169);
imagesc(X,Y,u_xy);
title('(e) Im(Ug)','FontSize',11,'fontweight','bold');
xlabel('x(m)');
ylabel('y(m)');
view(2);axis image;colormap jet;
caxis([-100,100]);
set(gca,'FontSize',11,'fontweight','bold');
subplot(3,3,6)   %------Absolute value------%
u_xy=reshape(full(abs(u(:))),169,169);
imagesc(X,Y,u_xy);
title('(f) Abs(Ug)','FontSize',11,'fontweight','bold');
xlabel('x(m)');
ylabel('y(m)');
view(2);axis image;colormap jet;
caxis([-100,100]);
set(gca,'FontSize',11,'fontweight','bold');
%----------------------------------------Optimal PML----------------------------------------------%
subplot(3,3,7)   %------Real part------%
u_new_xy=reshape(full(real(u_new(:))),169,169);
imagesc(X,Y,u_new_xy);
title('(g) Re(Uo)','FontSize',11,'fontweight','bold');
xlabel('x(m)');
ylabel('y(m)');
view(2);axis image;colormap jet;
caxis([-100,100]);
set(gca,'FontSize',11,'fontweight','bold');
subplot(3,3,8)   %------Imaginary part------%
u_new_xy=reshape(full(imag(u_new(:))),169,169);
imagesc(X,Y,u_new_xy);
title('(h) Im(Uo)','FontSize',11,'fontweight','bold');
xlabel('x(m)');
ylabel('y(m)');
view(2);axis image;colormap jet;
caxis([-100,100]);
set(gca,'FontSize',11,'fontweight','bold');
subplot9=subplot(3,3,9)   %------Absolute value------%
u_new_xy=reshape(full(abs(u_new(:))),169,169);
imagesc(X,Y,u_new_xy);
title('(i) Abs(Uo)','FontSize',11,'fontweight','bold');
xlabel('x(m)');
ylabel('y(m)');
view(2);axis image;colormap jet;
caxis([-100,100]);
colorbar('peer',subplot9,'Position',...
    [0.929853908619588 0.309551208285386 0.0201990623145487 0.356731872714326]);
set(gca,'FontSize',11,'fontweight','bold');
%% ------------Comparison maps of analytical solutions and two numerical solutions---------------%%
% -----------------------------comparation maps of wave value at different lines--------------------------------%
xx=X;               %x-coordinate
figure(4)
set(gcf,'Position',[113.0000 222.0000 1224 706.0000]);      %position of picture
% ---------------------------------x=(xl:ds1:xr),  y=0m-----------------------------------%
n_1=length(X);                 %Numbers of nodes of each line                
n_2=n_1*((n_1-1)*0.5)+1;       %the first node number at y=0m
Utrue1=Utrue(n_2:n_2+n_1-1);                %the line at y=0m
u1=u(n_2:n_2+n_1-1);
u2=u_new(n_2:n_2+n_1-1);
subplot(3,3,1)   %------Real part------%
plot(xx,real(Utrue1),'r-',xx,real(u1),'b:',xx,real(u2),'k--','LineWidth',1.8)
xlim([xl-d,xr+d]);
title('(a) Re(u) at y = 0 m','FontSize',12,'fontweight','bold');
xlabel('x(m)');
set(gca,'FontSize',12,'fontweight','bold');
subplot(3,3,4)   %------Imaginary part------%
plot(xx,imag(Utrue1),'r-',xx,imag(u1),'b:',xx,imag(u2),'k--','LineWidth',1.8)
xlim([xl-d,xr+d]);
title('(d) Im(u) at y = 0 m','FontSize',12,'fontweight','bold');
xlabel('x(m)');
set(gca,'FontSize',12,'fontweight','bold');
subplot(3,3,7)   %------Absolute value------%
plot(xx,abs(Utrue1),'r-',xx,abs(u1),'b:',xx,abs(u2),'k--','LineWidth',1.8)
xlim([xl-d,xr+d]);
title('(g) Abs(u) at y = 0 m','FontSize',12,'fontweight','bold');
xlabel('x(m)');
set(gca,'FontSize',12,'fontweight','bold');
% ---------------------------------------x=(xl:ds1:xr),  y=2m--------------------------------------------%
n_3=n_1*((n_1-1)*0.5+(xr/ds)*0.5)+1;   %the first node number at y=2m
Utrue1=Utrue(n_3:n_3+n_1-1);                        %the line at y=2m
u1=u(n_3:n_3+n_1-1);
u2=u_new(n_3:n_3+n_1-1);
subplot(3,3,2)   %------Real part------%
plot(xx,real(Utrue1),'r-',xx,real(u1),'b:',xx,real(u2),'k--','LineWidth',1.8)
xlim([xl-d,xr+d]);
ylim([-60,60]);
title('(b) Re(u) at y = 2 m','FontSize',12,'fontweight','bold');
xlabel('x(m)');
set(gca,'FontSize',12,'fontweight','bold');
subplot(3,3,5)   %------Imaginary part------%
plot(xx,imag(Utrue1),'r-',xx,imag(u1),'b:',xx,imag(u2),'k--','LineWidth',1.8)
xlim([xl-d,xr+d]);
ylim([-60,60]);
title('(e) Im(u) at y = 2 m','FontSize',12,'fontweight','bold');
xlabel('x(m)');
set(gca,'FontSize',12,'fontweight','bold');
subplot(3,3,8)   %------Absolute value------%
plot(xx,abs(Utrue1),'r-',xx,abs(u1),'b:',xx,abs(u2),'k--','LineWidth',1.8)
xlim([xl-d,xr+d]);
ylim([0,60]);
title('(h) Abs(u) at y = 2 m','FontSize',12,'fontweight','bold');
xlabel('x(m)');
set(gca,'FontSize',12,'fontweight','bold');
% --------------------------------------------x=(xl:ds1:xr),  y=4m-------------------------------------------%
n_4=n_1*((n_1-1)*0.5+xr/ds)+1;   %the first node number at y=4m
Utrue1=Utrue(n_4:n_4+n_1-1);                  %the line at y=4m
u1=u(n_4:n_4+n_1-1);
u2=u_new(n_4:n_4+n_1-1);
subplot(3,3,3)   %------Real part------%
plot(xx,real(Utrue1),'r-',xx,real(u1),'b:',xx,real(u2),'k--','LineWidth',1.8)
xlim([xl-d,xr+d]);
ylim([-60,60]);
title('(c) Re(u) at y = 4 m','FontSize',12,'fontweight','bold');
xlabel('x(m)');
set(gca,'FontSize',12,'fontweight','bold');
subplot(3,3,6)   %------Imaginary part------%
plot(xx,imag(Utrue1),'r-',xx,imag(u1),'b:',xx,imag(u2),'k--','LineWidth',1.8)
xlim([xl-d,xr+d]);
ylim([-60,60]);
title('(f) Im(u) at y = 4 m','FontSize',12,'fontweight','bold');
xlabel('x(m)');
set(gca,'FontSize',12,'fontweight','bold');
subplot(3,3,9)   %------Absolute value------%
plot(xx,abs(Utrue1),'r-',xx,abs(u1),'b:',xx,abs(u2),'k--','LineWidth',1.8)
xlim([xl-d,xr+d]);
ylim([0,60]);
legend1=legend('Analytical Solution','General PML','Optimal PML');
legend('boxoff');      
set(legend1,...
    'Position',[0.730936817422513 0.139754487685261 0.13725489886855 0.0849858333723404]);
set(legend1,'FontSize',12,'fontweight','bold')
title('(i) Abs(u) at y = 4 m','FontSize',12,'fontweight','bold');
xlabel('x(m)');
set(gca,'FontSize',12,'fontweight','bold');
%% --------Absolute error maps between analytical solutions and numerical solutions-------%%
figure(5)   
set(gcf,'Position',[236 186 983 624]);      %position of picture
%-----------------------analytical solution - numerical solution of general PML-------------------------------%
u_error4=(real(u)-real(Utrue));
u_error5=(imag(u)-imag(Utrue));
u_error6=(abs(u)-abs(Utrue));
u_error4(site)=0;u_error5(site)=0;u_error6(site)=0;  %the wave field value is extremely large, we set the error as 0
subplot(2,3,1)   %------Real part------%
u_error44=reshape(full(u_error4(:)),169,169);
imagesc(X,Y,u_error44);
title('(a) Re(Ug) - Re(Ut)','FontSize',12,'fontweight','bold');
xlabel('x(m)');
ylabel('y(m)');
view(2);axis image;colormap jet;
caxis([-5,5]);
xlim([-4,4]);ylim([-4,4]);
set(gca,'FontSize',12,'fontweight','bold');
subplot(2,3,2)   %------Imaginary part------%
u_error55=reshape(full(u_error5),169,169);
imagesc(X,Y,u_error55);
title('(b) Im(Ug) - Im(Ut)','FontSize',12,'fontweight','bold');
xlabel('x(m)');
ylabel('y(m)');
view(2);axis image;colormap jet;
caxis([-5,5]);
xlim([-4,4]);ylim([-4,4]);
set(gca,'FontSize',12,'fontweight','bold');
subplot(2,3,3)   %------Absolute value------%
u_error66=reshape(full(u_error6),169,169);
imagesc(X,Y,u_error66);
title('(c) Abs(Ug) - Abs(Ut)','FontSize',12,'fontweight','bold');
xlabel('x(m)');
ylabel('y(m)');
view(2);axis image;colormap jet;
caxis([-5,5]);
xlim([-4,4]);ylim([-4,4]);
set(gca,'FontSize',12,'fontweight','bold');
%-----------------------analytical solution - numerical solution of optimal PML-------------------------------%
u_error1=(real(u_new)-real(Utrue));  
u_error2=(imag(u_new)-imag(Utrue));  
u_error3=(abs(u_new)-abs(Utrue));  
u_error1(site)=0;u_error2(site)=0;u_error3(site)=0;  %the wave field value is extremely large, we set the error as 0
subplot(2,3,4)   %------Real part------%
u_error11=reshape(full(u_error1),169,169);
imagesc(X,Y,u_error11);
title('(d) Re(Uo) - Re(Ut)','FontSize',12,'fontweight','bold');
xlabel('x(m)');
ylabel('y(m)');
view(2);axis image;colormap jet;
caxis([-5,5]);
xlim([-4,4]);ylim([-4,4]);
set(gca,'FontSize',12,'fontweight','bold');
subplot(2,3,5)   %------Imaginary part------%
u_error22=reshape(full(u_error2),169,169);
imagesc(X,Y,u_error22);
title('(e) Im(Uo) - Im(Ut)','FontSize',12,'fontweight','bold');
xlabel('x(m)');
ylabel('y(m)');
view(2);axis image;colormap jet;
caxis([-5,5]);
xlim([-4,4]);ylim([-4,4]);
set(gca,'FontSize',12,'fontweight','bold');
subplot6=subplot(2,3,6);   %------Absolute value------%
u_error33=reshape(full(u_error3),169,169);
imagesc(X,Y,u_error33);
title('(f) Abs(Uo) - Abs(Ut)','FontSize',12,'fontweight','bold');
xlabel('x(m)');
ylabel('y(m)');
view(2);axis image;colormap jet;
caxis([-5,5]);
xlim([-4,4]);ylim([-4,4]);
colorbar('peer',subplot6,'Position',...
    [0.931278163673697 0.248397435897436 0.020358306188925 0.554074809371419]);
set(gca,'FontSize',12,'fontweight','bold');