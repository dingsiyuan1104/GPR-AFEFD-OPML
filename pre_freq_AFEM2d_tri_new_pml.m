function  forward_new=pre_freq_AFEM2d_tri_new_pml(mesh,pml_new,freq)
% Form structure for forward modeling with the optimal PML
% Useage
%    forward_new=pre_freq_AFEM2d_tri_new_pml(mesh,pml_new,freq)
% Input:
%         mesh  --------Number of nodes per triangular element
%      pml_new  --------Stracture of the optimal PML auxiliary parameters
%         freq  --------Frequency of forward modeling
% Output:
%  forward_new  --------Structure for forward modeling with the optimal PML
%
% -----------2018.10£¨Deshan Feng, Siyuan Ding, and Xun Wang£©----------- %
ep0=1/(36*pi)*1e-9;
miu=4*pi*1e-7;
ep = mesh.elem(:,4);
sig= mesh.elem(:,5);
omega = 2 * pi *freq;
k2 = omega^2.*ep0*miu*(ep-1i*sig/omega/ep0);
Sx =  1-(pml_new.sigx.*1i./sqrt(k2));
Sy =  1-(pml_new.sigy.*1i./sqrt(k2));
A = Sy./Sx;
B = Sx./Sy;
C = Sx.*Sy;
%% Element and Node
ND = size(mesh.node,1);
%% The area of the griangle element ( a & b )
a(:,1) = mesh.node(mesh.elem(:,2),2)-mesh.node(mesh.elem(:,3),2);    %a1=y2-y3
a(:,2) = mesh.node(mesh.elem(:,3),2)-mesh.node(mesh.elem(:,1),2);    %a2=y3-y1
a(:,3) = mesh.node(mesh.elem(:,1),2)-mesh.node(mesh.elem(:,2),2);    %a3=y1-y2
b(:,1) = mesh.node(mesh.elem(:,3),1)-mesh.node(mesh.elem(:,2),1);    %b1=x3-x2
b(:,2) = mesh.node(mesh.elem(:,1),1)-mesh.node(mesh.elem(:,3),1);    %b2=x1-x3
b(:,3) = mesh.node(mesh.elem(:,2),1)-mesh.node(mesh.elem(:,1),1);    %b3=x2-x1
area = 0.5*abs(a(:,1).*b(:,2)-a(:,2).*b(:,1));   %s=1/2*(a1*b2-a2*b1)
%% General synthetic system matrix
Me = [2,1,1;1,2,1;1,1,2]/12; %Element matrix
M = sparse(ND,ND);           %System matrices
K = sparse(ND,ND);           %System matrices
for i= 1:3
    for j = 1:3
        Mij = Me(i,j).*k2.*C.*area ;
        Kxij = a(:,i).*a(:,j)./(4*area).*A;
        Kyij = b(:,i).*b(:,j)./(4*area).*B;    
        Kij=Kxij+Kyij;
        
        M = M + sparse(mesh.elem(:,i),mesh.elem(:,j),Mij,ND,ND);
        K = K + sparse(mesh.elem(:,i),mesh.elem(:,j),Kij,ND,ND);
    end
end
K = M - K;         
forward_new.K=K;