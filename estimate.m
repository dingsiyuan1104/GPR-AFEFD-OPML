function eta = estimate(mesh)
% ESTIMATE computes ZZ error estimator on each element
% 
% USAGE
%    eta = estimate(mesh)
%
% INPUT 
%    mesh:  current mesh
%       u:  finite element solution
%
% OUTPUT
%    beta:  error estimator on each element
%

% L. Chen & C. Zhang 11-15-2006

[dudx,dudy] = ZZ(mesh);
[dudxx,dudxy] = sgrad(mesh,dudx);
[dudyx,dudyy] = sgrad(mesh,dudy);
eta = abs(dudxx) + abs(dudxy) + abs(dudyx) + abs(dudyy);

%--------------------------------------------------------------------------
% End of function ESTIMATE
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sub functions called by ESTIMATE
%--------------------------------------------------------------------------

function [dudx,dudy] = ZZ(mesh)
% ZZ is the Z-Z recovery error estimator  

N = size(mesh.node,1); 
NT = size(mesh.elem,1); 
activeNode = find(mesh.type~=0);   
nonActive = setdiff(1:N,activeNode);

[sdudx,sdudy,area] = sgrad(mesh,mesh.solu);
patcharea = sparse(mesh.elem(:,[1 2 3]),ones(NT,3),area*[1,1,1],N,1);

elem_0=mesh.elem(:,[1 2 3]);
dudx = accumarray(elem_0(:),[sdudx;sdudx;sdudx],[N 1]);      
dudy = accumarray(elem_0(:),[sdudy;sdudy;sdudy],[N 1]);      
dudx(activeNode) = dudx(activeNode)./patcharea(activeNode);
dudy(activeNode) = dudy(activeNode)./patcharea(activeNode);
dudx(nonActive) = 0; 
dudy(nonActive) = 0;
%--------------------------------------------------------------------------
% End of function ZZ
%--------------------------------------------------------------------------

function [sdudx,sdudy,area] = sgrad(mesh,u)
% SGRAD compute the scaled grad of w on the mesh
% it is elementwise
ve(:,:,1) = mesh.node(mesh.elem(:,3),:)-mesh.node(mesh.elem(:,2),:);
ve(:,:,2) = mesh.node(mesh.elem(:,1),:)-mesh.node(mesh.elem(:,3),:);
ve(:,:,3) = mesh.node(mesh.elem(:,2),:)-mesh.node(mesh.elem(:,1),:);
area = 0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));  
sdudx = -0.5*( u(mesh.elem(:,1)).*ve(:,2,1) ...
             + u(mesh.elem(:,2)).*ve(:,2,2) ...
             + u(mesh.elem(:,3)).*ve(:,2,3) );                
sdudy =  0.5*( u(mesh.elem(:,1)).*ve(:,1,1) ...
             + u(mesh.elem(:,2)).*ve(:,1,2) ...
             + u(mesh.elem(:,3)).*ve(:,1,3) );                
%--------------------------------------------------------------------------
% End of function SGRAD 
%--------------------------------------------------------------------------