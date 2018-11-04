%  2D Adaptive Finite Element Frequency Domain Simulation for GPR with PML 
% ---------------the optimal PML-------------unstructured mesh----------- %    
%
% 2D transverse magnetic(TM) mode waves
% d(dEz/dx)/dx+d(dEz/dy)/dy+k0^2*(epr-1i*sig/omega/eps0)*Ez=1i*omega*miu0*J;
%       k0^2 = omega^2*eps0*miu0;
%      d(A*du/dx)/dx+d(A*du/dy)/dy+C*k^2*u=0;
%      k^2 = omega^2*miu*eps0(epr-1i*sig/omega/eps0);
%      A = Sy/Sx ; B = Sx/Sy ; C = Sx*Sy
% optimal PML
%      Sx = 1 - 1i*sig_x/k;
%      sig_x = 1/(d-x);
%
% -----------2018.10£¨Deshan Feng, Siyuan Ding, and Xun Wang£©----------- %
close all;
clc;
clear;
addpath Solvers
load('outline_lining_fracture_model.mat')          %load the data
tic
%% Import COMSOL irregular triangulation mesh
[node, ND, elem, NT, elem_num, abc_node, abc_N, abc_num]= mphtxt_read_2d('mesh_optimal_PML_2_small.mphtxt');
%% ----------------------------Map of the model with initial triangulation------------------------------%
figure(1)
set(gcf,'Position',[131 435 562 390]);      %position of picture
i = 1:NT;
x = [node(elem(i,1),1) node(elem(i,2),1)  node(elem(i,3),1)]';
y = [node(elem(i,1),2) node(elem(i,2),2)  node(elem(i,3),2)]';
hs=patch('XData',x(:,i),'YData',y(:,i));
set(hs,'EdgeColor','w','LineWidth',0.5,'FaceColor','w');  % ---------PML-------%

ii = find(elem_num==7);   % ---------air layer--------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor','k','LineWidth',0.5,'FaceColor',[137, 207, 240]/255)
ii = find(elem_num==8);   % ---------air layer--------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor','k','LineWidth',0.5,'FaceColor',[137, 207, 240]/255)

ii = find(elem_num==12);  % ---------fracture--------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor',[137, 207, 240]/255,'LineWidth',0.5,'FaceColor','w')
ii = find(elem_num==13);  % ---------fracture--------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor',[137, 207, 240]/255,'LineWidth',0.5,'FaceColor','w')
ii = find(elem_num==14);  % ---------fracture--------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor',[137, 207, 240]/255,'LineWidth',0.5,'FaceColor','w')
ii = find(elem_num==15);  % ---------fracture--------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor',[137, 207, 240]/255,'LineWidth',0.5,'FaceColor','w')
ii = find(elem_num==16);  % ---------fracture--------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor',[137, 207, 240]/255,'LineWidth',0.5,'FaceColor','w')

ii = find(elem_num==9);   % ---------lining-------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor','k','LineWidth',0.5,'FaceColor',[193 255 193]/255)


ii = find(elem_num==10);  % ---------rock--------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor',[0.600000023841858 0.200000002980232 0],'LineWidth',0.5,'FaceColor',[231,179,37]/255)

axis ij
axis image
xlabel('Distance(m)')
ylabel('Depth(m)')

sitey = -0.01*ones(197,1);
sitex = 0.005:0.005:0.985;
plot(sitex,sitey,'y.','markersize',4)  %shot point
sitey = -0.01*ones(197,1);
recx = 0.015:0.005:0.995;
plot(recx,sitey,'rx','markersize',4)   %receivor

xlim([0,1]);ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                        
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...        
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                  
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'}); 
set(gca,'FontSize',12,'fontweight','bold');
%% Parameters of the medium
elem(:,4)=9;                      % relative permittivity---------lining-------%
elem(:,5)=0.00001;                % conductivity
elem(elem_num==10,4)=7;           % ---------rock--------%
elem(elem_num==10,5)=0.001;
elem(elem_num==7,4)=1;            % ---------air layer--------%
elem(elem_num==7,5)=0;
elem(elem_num==8,4)=1;            % ---------air layer--------%
elem(elem_num==8,5)=0;
elem(elem_num==12,4)=1;           % ---------fracture--------%
elem(elem_num==12,5)=0;
elem(elem_num==13,4)=1;           % ---------fracture--------%
elem(elem_num==13,5)=0;
elem(elem_num==14,4)=1;           % ---------fracture--------%
elem(elem_num==14,5)=0;
elem(elem_num==15,4)=1;           % ---------fracture--------%
elem(elem_num==15,5)=0;
elem(elem_num==16,4)=1;           % ---------fracture--------%
elem(elem_num==16,5)=0;
elem(elem_num==4,4)=7;                 % ---------PML around rock--------%
elem(elem_num==4,5)=0.001;
elem(elem_num==5,4)=7;                 % ---------PML around rock--------%
elem(elem_num==5,5)=0.001;
elem(elem_num==11,4)=7;                % ---------PML around rock--------%
elem(elem_num==11,5)=0.001;
elem(elem_num==20,4)=7;                % ---------PML around rock--------%
elem(elem_num==20,5)=0.001;
elem(elem_num==21,4)=7;                % ---------PML around rock--------%
elem(elem_num==21,5)=0.001;
elem(elem_num==1,4)=1;            % ---------PML around air layer--------%
elem(elem_num==1,5)=0;
elem(elem_num==2,4)=1;            % ---------PML around air layer--------%
elem(elem_num==2,5)=0;
elem(elem_num==6,4)=1;            % ---------PML around air layer--------%
elem(elem_num==6,5)=0;
elem(elem_num==17,4)=1;           % ---------PML around air layer--------%
elem(elem_num==17,5)=0;
elem(elem_num==18,4)=1;           % ---------PML around air layer--------%
elem(elem_num==18,5)=0;

ep0=1/(36*pi)*1e-9;
miu=4*pi*1e-7;
ep=elem(:,4);
sig=elem(:,5);
%% Load model parameters (area size and boundary partition)
xl = 0; xr = 1;         % The left and right boundaries of the calculation region in the X direction
yl = -0.1; yr = 0.5;    % The left and right boundaries of the calculation region in the Y direction
fprintf('the thickness of the optimal PML is:');
d =   0.03              % thickness of pml
%% Locations of the shot points and receiver points 
load('dsy_site_rec_xy_2.mat');        %x=0:0.005:1;   y=-0.01
sitexy = dsy_site_rec_xy_2(2:198,:);  %Locations of the shot points       x=0.005:0.005:0.985   
recxy = dsy_site_rec_xy_2(4:200,:);   %Locations of the receiver points   x=0.015:0.005:0.995
%% mesh parameters for AFEM
mesh0 = getmesh_AFEM_tri(node,elem,[],[]);             %node;elem;type;solu
%% Pre loading of excitation source
nsorc = size(sitexy,1);  %Number of the shot points
site = zeros(nsorc,1);
rec = zeros(nsorc,1);
for j = 1 : ND;
for i = 1 : nsorc;
    if abs(mesh0.node(j,:)-sitexy(i,:)) <= 0.0001
        site(i,1) = j;
    end
    if abs(mesh0.node(j,:)-recxy(i,:)) <= 0.0001
        rec(i,1) = j;
    end
end
end
%% Parameters of the antenna
freq =10e6:10e6:3000e6;       %----------Forward frequency-------------%
freq0 = 900e6;                % main frequency
%% ------------------------Ricker wavelet------------------------------%%
dt = 0.1*1e-9;    %  sampling period
fs = 1/dt;        %  sampling frequency
N = 500;          %  sampling numbers 
t = [0:N-1]*dt;  
x_t = (1-2*pi^2*(freq0*t-1).^2).*exp(-pi^2*(freq0*t-1).^2) ;    % Ricker wavelet in the time domain
%% Transformation from time domain to frequency domain 
Y = fft(x_t);          
P2 = Y(1:N/2+1)/N;  
nf = length(P2);
P1 = abs(P2);
f = fs*(0:(N/2))/N;           %----------Forward frequency-------------%
%% Inverse transformation from frequency domain to time domain
f_vol = zeros(N, 1);
f_vol(1:nf) = P2;
f_vol(nf+1:end) = conj(flip(P2(2:nf-1)));
data = ifft(N*f_vol);  
tt = [0:dt:dt*(length(data)-1)] ;  
data1 = real(data); 
%% ---------------------------plot------------------------------%%
%--------------------Waveform of the ricker wavelet in the time domain-----------------------% 
figure(2) 
subplot(311) 
plot(t/1e-9,x_t);  
title('Ricker');  xlabel('Tim(ns)');  ylabel('Amplitude'); 
%-------------------Spectrum of the ricker wavelet in the frequency domain-------------------% 
subplot(312)  
plot(f/1e6,P1);           
title('fft(Ricker)');  xlabel('Frequency(MHz)');  ylabel('Amplitude');
%-------Waveform of the inverse transformation from frequency domain to time domain----------%  
subplot(313)  
plot(tt/1e-9,data1);  
title('ifft(fft(Ricker))');  xlabel('Time(ns)');  ylabel('Amplitude');
%% Parameters for mesh refining
theta = 0.4; MaxIt = 2; RT = 3;
u_rec=zeros(nf,nsorc);% Record£¬which value in the first frequency is 0
[xq,yq]=meshgrid(-0.03:0.003:1.03,-0.13:0.003:0.53);% Interpolation grid spacing
u_rec_f = zeros(size(xq,1)*size(xq,2),nf);          % Record of one shot point£¬with interpolation mesh£¬the number of the sampling frequency is nf
for i=2:151          % main frequency, the loop starts from the second frequency
    i                % the sort number of the frequency
    omega = 2 * pi *f(i);
    for j=1:nsorc                       %ith shot point
%     for j=99                       %ith shot point
        %% Mesh Refinement
        mesh=mesh0;
        for iter = 1:MaxIt
            %% Step 1: Solve 
            ND = size(mesh.node,1);                       %Number of nodes
            RHS = zeros(ND,1);                            %Source vector
            RHS(site(j)) = 1i*omega*miu;                  
            pml_new = pml_new_freq_fem2d_tri(mesh.elem,mesh.node,xl,xr,yl,yr,d); % Calculate sigx and sigy
            forward_new=pre_freq_AFEM2d_tri_new_pml(mesh,pml_new,f(i));          % Generate the structure  
            K_new = forward_new.K;
            u_new=zeros(ND,1);
            u_new = opti_zmumps(K_new,RHS(:));
            mesh.solu=u_new(:);
            %% Step 2: Estimate   
            eta = estimate(mesh).^2;
            %% Step 3: Refine    
            for ii=1:RT
                [mesh,eta] = bisection(mesh,eta,theta);
            end
        end  % End of Mesh Refinement
        %% Step 4: Solve 
        ND = size(mesh.node,1);                       % Number of nodes
        RHS = zeros(ND,1);                            % Source vector
        RHS(site(j)) = 1i*omega*miu;                  
        pml_new = pml_new_freq_fem2d_tri(mesh.elem,mesh.node,xl,xr,yl,yr,d); % Calculate sigx and sigy
        forward_new=pre_freq_AFEM2d_tri_new_pml(mesh,pml_new,f(i));          % Generate the structure
        K_new = forward_new.K;
        u_new=zeros(ND,1);
        u_new = opti_zmumps(K_new,RHS(:));
        mesh.solu=u_new(:);
        u_rec(i,j)=u_new(rec(j))*P2(i);    %------------Record of jth shot point and ith frequency--------------%
        if j==99                           % Only save the record of jth shot point
            x1=mesh.node(:,1);
            y1=mesh.node(:,2);
            v =zeros(ND,1);
            v=(u_new)*P2(i);
            vq = zeros(size(xq,1),size(xq,2));
            vq(:,:)=griddata(x1,y1,v,xq,yq);                     % nodes of regular mesh after interpolation      
            u_rec_f(:,i) = reshape(vq,size(xq,1)*size(xq,2),1);  % jth shot point; ith frequency
            if i==21                              % Only save the adaptive mesh of 400MHz
                u_rec_f400(:)=u_new;
                mesh400=mesh;
            end
            if i==46                              % Only save the adaptive mesh of 900MHz
                u_rec_f900(:)=u_new;
                mesh900=mesh;
            end
            if i==76                              % Only save the adaptive mesh of 1500MHz
                u_rec_f1500(:)=u_new;
                mesh1500=mesh;
            end            
        end
    end  % End of loading the excitation source
end  % End the loop of frequency
toc
%% ------------------------------------Record ---------------------------------------%%
% -------------------Transform from time domain to frequency domain(ifft)-------------------------%
F_rec=zeros(N,nsorc);
F_rec(1:nf,:) = u_rec;                           %(N/2+1) values in front
F_rec(nf+1:end,:) = conj(flip(F_rec(2:nf-1,:))); %(N/2-1) values behind; turn over
data_rec = ifft(N*F_rec,[],1);                   %ifft 
tt = 0:dt:dt*(length(data_rec)-1) ; 
xx = 0.01:0.005:0.99;
data_rec = real(data_rec);

% single-channel waveform in the time domaim (99th shot point)
figure(3)
plot(tt*1e9, data_rec(:,99));
xlabel('Time(ns)');
ylabel('Amplitude');
xlim([0,15]);ylim([-2000,1500]);
set(gca,'FontSize',12,'fontweight','bold');

% profile in the time domain
% figure(5)
% imagesc(xx,tt*1e9,data_rec,[-200,200]);
% xlabel('Distance(m)');
% ylabel('Time(ns)');
% ylim([0,15]);
% set(gca,'FontSize',12,'fontweight','bold');
% colorbar;colormap jet;
%% ------------------ Wave field of jth shot point -------------------%
u_new_j=u_rec_f(:,:);    % N(interp) * N/2+1 
% -------------------Transform from time domain to frequency domain(ifft)-------------------------%
ff = f;%
F =  zeros(N,size(xq,1)*size(xq,2));
F1 = shiftdim(u_new_j,1);                % N/2+1 * N(interp)
F(1:nf,:) = F1;                          % (N/2+1) values in front
F(nf+1:end,:) = conj(flip(F1(2:nf-1,:)));% (N/2-1) values behind; turn over
data0 = ifft(N*F,[],1);                  % ifft
%% ----------Wave field snapshots of different time with optimal PML(shot point is located at 0.5m)------------%
figure(11)   
set(gcf,'Position',[503 457 1076 450]);      %position of picture
%------------------------------------ 1.2ns ---------------------------------------%
vq = full(real(data0(13,:)));
vq = reshape(vq,size(xq,1),size(xq,2));
% Fig.11(a)
subplot(2,3,1)
imagesc(xq(1,:),yq(:,1),vq);
hold on    %------plot the outline lining fracture and rock------ %
plot(frac1_x,frac1_y,'k--','LineWidth',0.5);
hold on
plot(frac2_x,frac2_y,'k--','LineWidth',0.5);
hold on
plot(rockx,rocky,'k--','LineWidth',0.5);
hold on
plot([0 1],[0 0],'k--','LineWidth',0.5);
title('(a) t=1.2ns','FontSize',12,'fontweight','bold');
xlabel('Distance(m)');
ylabel('Depth(m)');
axis image;colormap jet;
caxis([-200,200]);
xlim([0,1]);
ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                         
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...        
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                 
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'});
set(gca,'FontSize',12,'fontweight','bold');
%------------------------------------ 2.8ns ---------------------------------------%
vq = full(real(data0(29,:)));
vq = reshape(vq,size(xq,1),size(xq,2));
% Fig.11(b)
subplot(2,3,2)
imagesc(xq(1,:),yq(:,1),vq);
hold on    %------plot the outline lining fracture and rock------ %
plot(frac1_x,frac1_y,'k--','LineWidth',0.8);
hold on
plot(frac2_x,frac2_y,'k--','LineWidth',0.8);
hold on
plot(rockx,rocky,'k--','LineWidth',0.8);
hold on
plot([0 1],[0 0],'k--','LineWidth',0.8);
title('(b) t=2.8ns','FontSize',12,'fontweight','bold');
xlabel('Distance(m)');
ylabel('Depth(m)');
axis image;colormap jet;
caxis([-200,200]);
xlim([0,1]);
ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                        
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...       
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                 
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'}); 
set(gca,'FontSize',12,'fontweight','bold');
%------------------------------------ 3.8ns ---------------------------------------%
vq = full(real(data0(39,:)));
vq = reshape(vq,size(xq,1),size(xq,2));
% Fig.11(c)
subplot(2,3,3)
imagesc(xq(1,:),yq(:,1),vq);
hold on    %------plot the outline lining fracture and rock------ %
plot(frac1_x,frac1_y,'k--','LineWidth',0.8);
hold on
plot(frac2_x,frac2_y,'k--','LineWidth',0.8);
hold on
plot(rockx,rocky,'k--','LineWidth',0.8);
hold on
plot([0 1],[0 0],'k--','LineWidth',0.8);
title('(c) t=3.8ns','FontSize',12,'fontweight','bold');
xlabel('Distance(m)');
ylabel('Depth(m)');
axis image;colormap jet;
caxis([-200,200]);
xlim([0,1]);
ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                     
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...        
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'});
set(gca,'FontSize',12,'fontweight','bold');
%------------------------------------ 4.9ns ---------------------------------------%
vq = full(real(data0(50,:)));
vq = reshape(vq,size(xq,1),size(xq,2));
% Fig.11(d)
subplot(2,3,4)
imagesc(xq(1,:),yq(:,1),vq);
hold on    %------plot the outline lining fracture and rock------ %
plot(frac1_x,frac1_y,'k--','LineWidth',0.8);
hold on
plot(frac2_x,frac2_y,'k--','LineWidth',0.8);
hold on
plot(rockx,rocky,'k--','LineWidth',0.8);
hold on
plot([0 1],[0 0],'k--','LineWidth',0.8);
title('(d) t=4.9ns','FontSize',12,'fontweight','bold');
xlabel('Distance(m)');
ylabel('Depth(m)');
axis image;colormap jet;
caxis([-200,200]);
xlim([0,1]);
ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                       
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...      
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                 
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'}); 
set(gca,'FontSize',12,'fontweight','bold');
%------------------------------------ 6.3ns ---------------------------------------%
vq = full(real(data0(64,:)));
vq = reshape(vq,size(xq,1),size(xq,2));
% Fig.11(e)
subplot(2,3,5)
imagesc(xq(1,:),yq(:,1),vq);
hold on    %------plot the outline lining fracture and rock------ %
plot(frac1_x,frac1_y,'k--','LineWidth',0.8);
hold on
plot(frac2_x,frac2_y,'k--','LineWidth',0.8);
hold on
plot(rockx,rocky,'k--','LineWidth',0.8);
hold on
plot([0 1],[0 0],'k--','LineWidth',0.8);
title('(e) t=6.3ns','FontSize',12,'fontweight','bold');
xlabel('Distance(m)');
ylabel('Depth(m)');
axis image;colormap jet;
caxis([-200,200]);
xlim([0,1]);
ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                   
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...      
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                  
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'});
set(gca,'FontSize',12,'fontweight','bold');
%------------------------------------ 7.5ns ---------------------------------------%
vq = full(real(data0(76,:)));
vq = reshape(vq,size(xq,1),size(xq,2));
% Fig.11(f)
subplot6=subplot(2,3,6)
imagesc(xq(1,:),yq(:,1),vq);
hold on    %------plot the outline lining fracture and rock------ %
plot(frac1_x,frac1_y,'k--','LineWidth',0.8);
hold on
plot(frac2_x,frac2_y,'k--','LineWidth',0.8);
hold on
plot(rockx,rocky,'k--','LineWidth',0.8);
hold on
plot([0 1],[0 0],'k--','LineWidth',0.8);
title('(f) t=7.5ns','FontSize',12,'fontweight','bold');
xlabel('Distance(m)');
ylabel('Depth(m)');
axis image;colormap jet;
caxis([-200,200]);
xlim([0,1]);
ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                         
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...      
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'}); 
set(gca,'FontSize',12,'fontweight','bold');
colorbar('peer',subplot6,'Position',...
    [0.930083608932516 0.364444444444444 0.018626373626374 0.360575798811394]);
%% --Profiles of different frequency with the optimal PML(shot point is located at 0.5m(99th shot point))-- %%
figure(10)   
set(gcf,'Position',[76 137 1314 825]);      %position of picture
%-----------------------------------frequency of 400MHz-----------------------------------------%   
u_new1=u_rec_f400;
[xq,yq]=meshgrid(-0.03:0.0015:1.03,-0.13:0.0015:0.53);
x1=mesh400.node(:,1);
y1=mesh400.node(:,2);
v =full(u_new1);
vq = griddata(x1,y1,v,xq,yq);
% Fig.10(a)
subplot(3,3,1)                 %---------Real part----------%
imagesc(xq(1,:),yq(:,1),real(vq));

hold on    %------plot the outline lining fracture and rock------ %
plot(frac1_x,frac1_y,'k--','LineWidth',0.8);
hold on
plot(frac2_x,frac2_y,'k--','LineWidth',0.8);
hold on
plot(rockx,rocky,'k--','LineWidth',0.8);
hold on
plot([0 1],[0 0],'k--','LineWidth',0.8);

title('(a) Re(u) of 400MHz','FontSize',12,'fontweight','bold');
xlabel('Distance(m)');
ylabel('Depth(m)');
axis image;colormap jet;
caxis([-800,800]);
xlim([0,1]);
ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                        
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...       
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                 
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'}); 
set(gca,'FontSize',12,'fontweight','bold');
% Fig.10(b)
subplot(3,3,2)                 %---------Imaginary part----------%
imagesc(xq(1,:),yq(:,1),imag(vq));

hold on    %------plot the outline lining fracture and rock------ %
plot(frac1_x,frac1_y,'k--','LineWidth',0.8);
hold on
plot(frac2_x,frac2_y,'k--','LineWidth',0.8);
hold on
plot(rockx,rocky,'k--','LineWidth',0.8);
hold on
plot([0 1],[0 0],'k--','LineWidth',0.8);

title('(b) Im(u) of 400MHz','FontSize',12,'fontweight','bold');
xlabel('Distance(m)');
ylabel('Depth(m)');
axis image;colormap jet;
caxis([-800,800]);
xlim([0,1]);
ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                        
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...       
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                 
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'}); 
set(gca,'FontSize',12,'fontweight','bold');
% Fig.10(c)
subplot(3,3,3)                 %---------Absolute value----------%
imagesc(xq(1,:),yq(:,1),abs(vq));
hold on    %------plot the outline lining fracture and rock------ %
plot(frac1_x,frac1_y,'k--','LineWidth',0.8);
hold on
plot(frac2_x,frac2_y,'k--','LineWidth',0.8);
hold on
plot(rockx,rocky,'k--','LineWidth',0.8);
hold on
plot([0 1],[0 0],'k--','LineWidth',0.8);
title('(c) Abs(u) of 400MHz','FontSize',12,'fontweight','bold');
xlabel('Distance(m)');
ylabel('Depth(m)');
axis image;colormap jet;
caxis([-800,800]);
xlim([0,1]);
ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                         
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...        
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                  
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'}); 
set(gca,'FontSize',12,'fontweight','bold');
%-----------------------------------frequency of 900MHz-----------------------------------------%  
u_new1=u_rec_f900;   
x1=mesh900.node(:,1);
y1=mesh900.node(:,2);
v =full(u_new1);
vq = griddata(x1,y1,v,xq,yq);
% Fig.10(d)
subplot(3,3,4)                 %---------Real part----------%
imagesc(xq(1,:),yq(:,1),real(vq));
hold on    %------plot the outline lining fracture and rock------ %
plot(frac1_x,frac1_y,'k--','LineWidth',0.8);
hold on
plot(frac2_x,frac2_y,'k--','LineWidth',0.8);
hold on
plot(rockx,rocky,'k--','LineWidth',0.8);
hold on
plot([0 1],[0 0],'k--','LineWidth',0.8);
title('(d) Re(u) of 900MHz','FontSize',12,'fontweight','bold');
xlabel('Distance(m)');
ylabel('Depth(m)');
axis image;colormap jet;
caxis([-800,800]);
xlim([0,1]);
ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                        
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...       
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                  
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'}); 
set(gca,'FontSize',12,'fontweight','bold');
% Fig.10(e)
subplot(3,3,5)                 %---------Imaginary part----------%
imagesc(xq(1,:),yq(:,1),imag(vq));
hold on    %------plot the outline lining fracture and rock------ %
plot(frac1_x,frac1_y,'k--','LineWidth',0.8);
hold on
plot(frac2_x,frac2_y,'k--','LineWidth',0.8);
hold on
plot(rockx,rocky,'k--','LineWidth',0.8);
hold on
plot([0 1],[0 0],'k--','LineWidth',0.8);
title('(e) Im(u) of 900MHz','FontSize',12,'fontweight','bold');
xlabel('Distance(m)');
ylabel('Depth(m)');
axis image;colormap jet;
caxis([-800,800]);
xlim([0,1]);
ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                        
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...       
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                 
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'});
set(gca,'FontSize',12,'fontweight','bold');
% Fig.10(f)
subplot(3,3,6)                 %---------Absolute value----------%
imagesc(xq(1,:),yq(:,1),abs(vq));
hold on    %------plot the outline lining fracture and rock------ %
plot(frac1_x,frac1_y,'k--','LineWidth',0.8);
hold on
plot(frac2_x,frac2_y,'k--','LineWidth',0.8);
hold on
plot(rockx,rocky,'k--','LineWidth',0.8);
hold on
plot([0 1],[0 0],'k--','LineWidth',0.8);
title('(f) Abs(u) of 900MHz','FontSize',12,'fontweight','bold');
xlabel('Distance(m)');
ylabel('Depth(m)');
axis image;colormap jet;
caxis([-800,800]);
xlim([0,1]);
ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                        
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...       
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                 
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'}); 
set(gca,'FontSize',12,'fontweight','bold');
%-----------------------------------frequency of 1500MHz-----------------------------------------%
u_new1=u_rec_f1500; 
x1=mesh1500.node(:,1);
y1=mesh1500.node(:,2);
v =full(u_new1);
vq = griddata(x1,y1,v,xq,yq);
% Fig.10(g)
subplot(3,3,7)                 %---------Real part----------%
imagesc(xq(1,:),yq(:,1),real(vq));
hold on    %------plot the outline lining fracture and rock------ %
plot(frac1_x,frac1_y,'k--','LineWidth',0.8);
hold on
plot(frac2_x,frac2_y,'k--','LineWidth',0.8);
hold on
plot(rockx,rocky,'k--','LineWidth',0.8);
hold on
plot([0 1],[0 0],'k--','LineWidth',0.8);
title('(g) Re(u) of 1500MHz','FontSize',12,'fontweight','bold');
xlabel('Distance(m)');
ylabel('Depth(m)');
axis image;colormap jet;
caxis([-800,800]);
xlim([0,1]);
ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                        
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...        
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                  
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'}); 
set(gca,'FontSize',12,'fontweight','bold');
% Fig.10(h)
subplot(3,3,8)                 %---------Imaginary part----------%
imagesc(xq(1,:),yq(:,1),imag(vq));
hold on    %------plot the outline lining fracture and rock------ %
plot(frac1_x,frac1_y,'k--','LineWidth',0.8);
hold on
plot(frac2_x,frac2_y,'k--','LineWidth',0.8);
hold on
plot(rockx,rocky,'k--','LineWidth',0.8);
hold on
plot([0 1],[0 0],'k--','LineWidth',0.8);
title('(h) Im(u) of 1500MHz','FontSize',12,'fontweight','bold');
xlabel('Distance(m)');
ylabel('Depth(m)');
axis image;colormap jet;
caxis([-800,800]);
xlim([0,1]);
ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                         
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...        
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                  
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'}); 
set(gca,'FontSize',12,'fontweight','bold');
% Fig.10(i)
subplot9=subplot(3,3,9)                 %---------Absolute value----------%
imagesc(xq(1,:),yq(:,1),abs(vq));
hold on    %------plot the outline lining fracture and rock------ %
plot(frac1_x,frac1_y,'k--','LineWidth',0.8);
hold on
plot(frac2_x,frac2_y,'k--','LineWidth',0.8);
hold on
plot(rockx,rocky,'k--','LineWidth',0.8);
hold on
plot([0 1],[0 0],'k--','LineWidth',0.8);
title('(i) Abs(u) at 1500MHz','FontSize',12,'fontweight','bold');
xlabel('Distance(m)');
ylabel('Depth(m)');
axis image;colormap jet;
caxis([-800,800]);
xlim([0,1]);
ylim([-0.1,0.5]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],...                         
        'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...       
        'YTick',[-0.1 0 0.1 0.2 0.3 0.4 0.5],...                  
        'YTickLabel',{'-0.1','0','0.1','0.2','0.3','0.4','0.5'}); 
set(gca,'FontSize',12,'fontweight','bold');
colorbar('peer',subplot9,'Position',...
    [0.929056047078958 0.315716272600834 0.0224336283185844 0.36954102955748]);
%% ------------------------------Adaptive meshes(AFEFM in MATLAB)--------------------------------------%
figure(9)
set(gcf,'Position',[163 119 968 809]);      %position of picture
%% -----------------------------------grequency of 400MHz-----------------------------------------%
subplot(1,3,1)                
NT = length(mesh400.elem(:,1));
i = 1:NT;
x = [mesh400.node(mesh400.elem(i,1),1)  mesh400.node(mesh400.elem(i,2),1)  mesh400.node(mesh400.elem(i,3),1)]';
y = [mesh400.node(mesh400.elem(i,1),2)  mesh400.node(mesh400.elem(i,2),2)  mesh400.node(mesh400.elem(i,3),2)]';
hs=patch('XData',x(:,i),'YData',y(:,i));
set(hs,'EdgeColor','w','LineWidth',0.5,'FaceColor','w');  % ---------PML-------%

ii = find(mesh400.elem(:,4)==1);   % ---------air layer and fracture--------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor','k','LineWidth',0.5,'FaceColor',[137, 207, 240]/255)

xx=(mesh400.node(mesh400.elem(:,1),1)+mesh400.node(mesh400.elem(:,2),1)+mesh400.node(mesh400.elem(:,3),1))/3;
yy=(mesh400.node(mesh400.elem(:,1),2)+mesh400.node(mesh400.elem(:,2),2)+mesh400.node(mesh400.elem(:,3),2))/3;
[in_1,~] = inpolygon(xx,yy,frac1_x,frac1_y);        %fracture_1
[in_2,~] = inpolygon(xx,yy,frac2_x,frac2_y);        %fracture_2
hold on
hs=patch('XData',x(:,in_1),'YData',y(:,in_1));% ---------fracture_1--------%
set(hs,'EdgeColor',[137, 207, 240]/255,'LineWidth',0.5,'FaceColor','w')
hold on
hs=patch('XData',x(:,in_2),'YData',y(:,in_2));% ---------fracture_2--------%
set(hs,'EdgeColor',[137, 207, 240]/255,'LineWidth',0.5,'FaceColor','w')

ii = find(mesh400.elem(:,4)==9);   % ---------lining-------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor','k','LineWidth',0.5,'FaceColor',[193 255 193]/255)

ii = find(mesh400.elem(:,4)==7);  % ---------rock--------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor',[0.600000023841858 0.200000002980232 0],'LineWidth',0.5,'FaceColor',[231,179,37]/255)

axis ij
axis image
% axis off
xlabel('Distance(m)')
ylabel('Depth(m)')

sitey = -0.01;
sitex = 0.5;
plot(sitex,sitey,'y.','markersize',8)  %shot point
xlim([0.4,0.6]);ylim([-0.05,0.45]);


set(gca,'XTick',[0.4  0.5 0.6],...                         
        'XTickLabel',{'0.4','0.5','0.6'},...       
        'YTick',[0 0.1 0.2 0.3 0.4],...                 
        'YTickLabel',{'0','0.1','0.2','0.3','0.4'}); 
title('(a) 400MHz','FontSize',12,'fontweight','bold');
set(gca,'FontSize',12,'fontweight','bold');
%% -----------------------------------frequency of 900MHz-----------------------------------------%
subplot(1,3,2)                
NT = length(mesh900.elem(:,1));
i = 1:NT;
x = [mesh900.node(mesh900.elem(i,1),1)  mesh900.node(mesh900.elem(i,2),1)  mesh900.node(mesh900.elem(i,3),1)]';
y = [mesh900.node(mesh900.elem(i,1),2)  mesh900.node(mesh900.elem(i,2),2)  mesh900.node(mesh900.elem(i,3),2)]';
hs=patch('XData',x(:,i),'YData',y(:,i));
set(hs,'EdgeColor','w','LineWidth',0.5,'FaceColor','w');  % ---------PML-------%

ii = find(mesh900.elem(:,4)==1);   % ---------air layer and fracture--------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor','k','LineWidth',0.5,'FaceColor',[137, 207, 240]/255)

xx=(mesh900.node(mesh900.elem(:,1),1)+mesh900.node(mesh900.elem(:,2),1)+mesh900.node(mesh900.elem(:,3),1))/3;
yy=(mesh900.node(mesh900.elem(:,1),2)+mesh900.node(mesh900.elem(:,2),2)+mesh900.node(mesh900.elem(:,3),2))/3;
[in_1,~] = inpolygon(xx,yy,frac1_x,frac1_y);        %fracture_1
[in_2,~] = inpolygon(xx,yy,frac2_x,frac2_y);        %fracture_2
hold on
hs=patch('XData',x(:,in_1),'YData',y(:,in_1));% ---------fracture_1--------%
set(hs,'EdgeColor',[137, 207, 240]/255,'LineWidth',0.5,'FaceColor','w')
hold on
hs=patch('XData',x(:,in_2),'YData',y(:,in_2));% ---------fracture_2--------%
set(hs,'EdgeColor',[137, 207, 240]/255,'LineWidth',0.5,'FaceColor','w')

ii = find(mesh900.elem(:,4)==9);   % ---------lining-------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor','k','LineWidth',0.5,'FaceColor',[193 255 193]/255)


ii = find(mesh900.elem(:,4)==7);  % ---------rock--------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor',[0.600000023841858 0.200000002980232 0],'LineWidth',0.5,'FaceColor',[231,179,37]/255)

axis ij
axis image
% axis off
xlabel('Distance(m)')
ylabel('Depth(m)')

sitey = -0.01;
sitex = 0.5;
plot(sitex,sitey,'y.','markersize',8)  %shot point
xlim([0.4,0.6]);ylim([-0.05,0.45]);


set(gca,'XTick',[0.4  0.5 0.6],...                         
        'XTickLabel',{'0.4','0.5','0.6'},...       
        'YTick',[0 0.1 0.2 0.3 0.4],...                 
        'YTickLabel',{'0','0.1','0.2','0.3','0.4'}); 
title('(b) 900MHz','FontSize',12,'fontweight','bold');
set(gca,'FontSize',12,'fontweight','bold');
%% -----------------------------------frequency of 1500MHz-----------------------------------------%
subplot(1,3,3)                 
NT = length(mesh1500.elem(:,1));
i = 1:NT;
x = [mesh1500.node(mesh1500.elem(i,1),1)  mesh1500.node(mesh1500.elem(i,2),1)  mesh1500.node(mesh1500.elem(i,3),1)]';
y = [mesh1500.node(mesh1500.elem(i,1),2)  mesh1500.node(mesh1500.elem(i,2),2)  mesh1500.node(mesh1500.elem(i,3),2)]';
hs=patch('XData',x(:,i),'YData',y(:,i));
set(hs,'EdgeColor','w','LineWidth',0.5,'FaceColor','w');  % ---------PML-------%

ii = find(mesh1500.elem(:,4)==1);   % ---------air layer and fracture--------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor','k','LineWidth',0.5,'FaceColor',[137, 207, 240]/255)

xx=(mesh1500.node(mesh1500.elem(:,1),1)+mesh1500.node(mesh1500.elem(:,2),1)+mesh1500.node(mesh1500.elem(:,3),1))/3;
yy=(mesh1500.node(mesh1500.elem(:,1),2)+mesh1500.node(mesh1500.elem(:,2),2)+mesh1500.node(mesh1500.elem(:,3),2))/3;
[in_1,~] = inpolygon(xx,yy,frac1_x,frac1_y);        %fracture_1
[in_2,~] = inpolygon(xx,yy,frac2_x,frac2_y);        %fracture_2
hold on
hs=patch('XData',x(:,in_1),'YData',y(:,in_1));% ---------fracture_1--------%
set(hs,'EdgeColor',[137, 207, 240]/255,'LineWidth',0.5,'FaceColor','w')
hold on
hs=patch('XData',x(:,in_2),'YData',y(:,in_2));% ---------fracture_2--------%
set(hs,'EdgeColor',[137, 207, 240]/255,'LineWidth',0.5,'FaceColor','w')

ii = find(mesh1500.elem(:,4)==9);   % ---------lining-------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor','k','LineWidth',0.5,'FaceColor',[193 255 193]/255)


ii = find(mesh1500.elem(:,4)==7);  % ---------rock--------%
hold on
hs=patch('XData',x(:,ii),'YData',y(:,ii));
set(hs,'EdgeColor',[0.600000023841858 0.200000002980232 0],'LineWidth',0.5,'FaceColor',[231,179,37]/255)

axis ij
axis image

xlabel('Distance(m)')
ylabel('Depth(m)')

sitey = -0.01;
sitex = 0.5;
plot(sitex,sitey,'y.','markersize',8)  %shot point
xlim([0.4,0.6]);ylim([-0.05,0.45]);


set(gca,'XTick',[0.4  0.5 0.6],...                         
        'XTickLabel',{'0.4','0.5','0.6'},...       
        'YTick',[0 0.1 0.2 0.3 0.4],...                 
        'YTickLabel',{'0','0.1','0.2','0.3','0.4'}); 
title('(c) 1500MHz','FontSize',12,'fontweight','bold');
set(gca,'FontSize',12,'fontweight','bold');