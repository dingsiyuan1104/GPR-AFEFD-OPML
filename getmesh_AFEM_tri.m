function mesh=getmesh_AFEM_tri(node,elem,Dirichlet,Neumann)
% Generate the grid data for forward modeling (triangular mesh)
% Useage
%    mesh = getmesh_AFEM_tri(node,elem,Dirichlet,Neumann)
% Input:
%         node  --------Mesh point coordinates
%         elem  --------Number of nodes per triangular element
%     Dirchlet  --------The first boundary condition (natural boundary condition)
%      Neumann  --------The third kind of boundary conditions (ABC)
% Output:
%         mesh  --------Mesh structure
%
% -----------2018.10£¨Deshan Feng, Siyuan Ding, and Xun Wang£©----------- %
if nargin <= 4, Dirichlet=[]; end
if nargin <= 3,  Neumann=[]; end
N=size(node,1);    
type=ones(N,1);
solu=zeros(N,1);

mesh.node=node;           
mesh.elem=elem;           
mesh.Dirichlet=Dirichlet; 
mesh.Neumann=Neumann;     
mesh.type=type;           % Node type
mesh.solu=solu;           % Solution 