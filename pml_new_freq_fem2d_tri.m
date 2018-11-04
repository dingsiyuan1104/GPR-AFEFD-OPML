function pml_new=pml_new_freq_fem2d_tri(elem,node,xl,xr,yl,yr,d)
% Calculate the optimal PML auxiliary parameters in the frequency domain forward (sigx,sigy)
% Useage
%    pml_new=pml_new_freq_fem2d_tri(elem,node,xl,xr,yl,yr,d)
% Input:
%     elem  --------Number of nodes per triangular element
%     node  --------Number of mesh points
%       xl  --------The left boundary of the calculation region in the X direction
%       xr  --------The right boundary of the calculation region in the X direction
%       yl  --------The left boundary of the calculation region in the Y direction
%       yr  --------The right boundary of the calculation region in the Y direction
%        d  --------Thickness of pml
% Output:
%  pml_new  --------Stracture of the optimal PML auxiliary parameters
%
% -----------2018.10（Deshan Feng, Siyuan Ding, and Xun Wang）----------- %
% ep0=1/(36*pi)*1e-9;
% miu=4*pi*1e-7; 
ep = elem(:,4);
% the cordinates of element in the X & Y directions
elem_x=(node(elem(:,1),1)+node(elem(:,2),1)+node(elem(:,3),1))/3;%三个顶点x坐标平均值
elem_y=(node(elem(:,1),2)+node(elem(:,2),2)+node(elem(:,3),2))/3;%三个顶点y坐标平均值
% sigx
sigx = zeros(size(ep));
%  the left PML(in the X directions)
pmlx_l = find(elem_x<=xl);
sigx(pmlx_l)= 1./(d- (xl-elem_x(pmlx_l)));
%  the right PML(in the X directions)
pmlx_r = find(elem_x>=xr);
sigx(pmlx_r)= 1./(d- (elem_x(pmlx_r)-xr));
% sigy
sigy = zeros(size(ep));
%  the left PML(in the Y directions)
pmly_l = find(elem_y<=yl);
sigy(pmly_l)= 1./(d- (yl-elem_y(pmly_l)));
%  the right PML(in the Y directions)
pmly_r = find(elem_y>=yr);
sigy(pmly_r)= 1./(d- (elem_y(pmly_r)-yr));
% nodes of PML
node_px = find(node(:,1)<=xl);
node_px = [node_px;  find(node(:,1)>=xr)];
node_py = find(node(:,2)<=yl);
node_py = [node_py;  find(node(:,2)>=yr)];
% Form structure of the optimal PML
pml_new.sigx = sigx;
pml_new.sigy = sigy;
pml_new.node_px=node_px;
pml_new.node_py=node_py;