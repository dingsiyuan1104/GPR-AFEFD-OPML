function pml=pml_freq_fem2d_tri(elem,node,m,R,xl,xr,yl,yr,d)
% Calculate the general PML auxiliary parameters in the frequency domain forward (sigx,sigy)
% Useage
%    pml=pml_freq_fem2d_tri(elem,node,m,R,xl,xr,yl,yr,d)
% Input:
%     elem  --------Number of nodes per triangular element
%     node  --------Number of mesh points
%        m  --------Degree of polynomial grading (usually between 3 and 4)
%        R  --------Theoretical reflection coefficient of the target with normal incidence
%       xl  --------The left boundary of the calculation region in the X direction
%       xr  --------The right boundary of the calculation region in the X direction
%       yl  --------The left boundary of the calculation region in the Y direction
%       yr  --------The right boundary of the calculation region in the Y direction
%        d  --------Thickness of pml
% Output:
%      pml  --------Stracture of the general PML auxiliary parameters
%
% -----------2018.10£¨Deshan Feng, Siyuan Ding, and Xun Wang£©----------- %
% ep0=1/(36*pi)*1e-9;
% miu=4*pi*1e-7; 
ep = elem(:,4);
eta0=120*pi;
sigmax =-(m+1)*log(R)./(2*d*eta0*sqrt(ep));
% the cordinates of element in the X & Y directions
elem_x=(node(elem(:,1),1)+node(elem(:,2),1)+node(elem(:,3),1))/3;
elem_y=(node(elem(:,1),2)+node(elem(:,2),2)+node(elem(:,3),2))/3;
% sigx
sigx = zeros(size(ep));
%  the left PML(in the X directions)
pmlx_l = find(elem_x<=xl);
sigx(pmlx_l)= sigmax(pmlx_l).* (xl-elem_x(pmlx_l)).^m/d^m;
%  the right PML(in the X directions)
pmlx_r = find(elem_x>=xr);
sigx(pmlx_r)= sigmax(pmlx_r).* (elem_x(pmlx_r)-xr).^m/d^m;
% sigy
sigy = zeros(size(ep));
%  the left PML(in the Y directions)
pmly_l = find(elem_y<=yl);
sigy(pmly_l)= sigmax(pmly_l).* (yl-elem_y(pmly_l)).^m/d^m;
%  the right PML(in the Y directions)
pmly_r = find(elem_y>=yr);
sigy(pmly_r)= sigmax(pmly_r).* (elem_y(pmly_r)-yr).^m/d^m;
% nodes of PML
node_px = find(node(:,1)<=xl);
node_px = [node_px;  find(node(:,1)>=xr)];
node_py = find(node(:,2)<=yl);
node_py = [node_py;  find(node(:,2)>=yr)];
% Form structure of the general PML
pml.sigx = sigx;
pml.sigy = sigy;
pml.node_px=node_px;
pml.node_py=node_py;