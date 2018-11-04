function [nodes, nodeSize, triangle, triangleSize, trianglebelong, surface, surfaceSize, surfacebelong]=mphtxt_read_2d(MeshFile)
% *.mphtxt file read
% Read COMSOL triangulation result
% Useage
%    [node, ND, elem, NT, elem_num, abc_node, abc_N, abc_num]= mphtxt_read_2d('name_of_file.mphtxt')
% Input:
% name_of_file.mphtxt  --------file of COMSOL triangulation result
% Output:
%               nodes  --------Mesh point coordinates
%            nodeSize  --------Number of mesh points
%            triangle  --------Number of nodes per triangular element
%        triangleSize  --------Number of triangular elements
%      trianglebelong  --------Number of geometric entity indices (refers to triangular element)
%             surface  --------Number of nodes per edge element
%         surfaceSize  --------Number of edge elements
%       surfacebelong  --------Number of geometric entity indices (refers to edge element)
%
% -----------2018.10£¨Deshan Feng, Siyuan Ding, and Xun Wang£©----------- %
fp = fopen(MeshFile);

flag=0;
count=0;
while (flag==0)&&(count<1000000);
   temp=fgetl(fp);         
   temp1='2 # sdim';
   flag=strcmp(temp, temp1);
   if flag==1
       temp2='# number of mesh points';
       nodeSize=fscanf(fp,['%d' temp2]);
   end
   count=count+1; 
end


flag=0;
count=0;
while (flag==0)&&(count<1000000);
   temp=fgetl(fp);
   temp1='# Mesh point coordinates';
   flag=strcmp(temp, temp1);
   if flag==1
       nodes=fscanf(fp,'%f %f',[2,nodeSize]);
   end
   count=count+1; 
end
nodes=nodes';          


flag=0;
count=0;
while (flag==0)&&(count<1000000);
   temp=fgetl(fp);
   temp1='2 # number of nodes per element';
   flag=strcmp(temp, temp1);
   if flag==1
       temp2='# number of elements';
       surfaceSize=fscanf(fp,['%d' temp2]);
       fgetl(fp);
       fgetl(fp);
       surface=fscanf(fp,'%f %f',[2,surfaceSize]);
   end
   count=count+1; 
end
surface=surface'+1;



flag=0;
count=0;
while (flag==0)&&(count<1000000);
   temp=fgetl(fp);
   temp1=[num2str(surfaceSize) ' # number of geometric entity indices'];
   flag=strcmp(temp, temp1);
   if flag==1
       fgetl(fp);
       surfacebelong=fscanf(fp,'%f',[1,surfaceSize]);
   end
   count=count+1; 
end
surfacebelong=surfacebelong'+1;


flag=0;
count=0;
while (flag==0)&&(count<1000000);
   temp=fgetl(fp);
   temp1='3 # number of nodes per element';
   flag=strcmp(temp, temp1);
   if flag==1
       temp2='# number of elements';
       triangleSize=fscanf(fp,['%d' temp2]);
       fgetl(fp);
       fgetl(fp);
       triangle=fscanf(fp,'%f %f %f',[3,triangleSize]);
   end
   count=count+1; 
end
triangle=triangle'+1;


flag=0;
count=0;
while (flag==0)&&(count<1000000);
   temp=fgetl(fp);
   temp1=[num2str(triangleSize) ' # number of geometric entity indices'];
   flag=strcmp(temp, temp1);
   if flag==1
       fgetl(fp);
       trianglebelong=fscanf(fp,'%f',[1,triangleSize]);
   end
   count=count+1; 
end
trianglebelong=trianglebelong';

fclose(fp);