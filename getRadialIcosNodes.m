function [X, tri] = getRadialIcosNodes(m,n)
%%Returns Icosahedral nodes. See Caspar, Klug. Projection is radial
% Will return total points T.
%
%   T = 10(m^2+mn+n^2) + 2
%
%
% [X,tri] returns a triangulation of the nodes.
%
% Author: T. Michaels
% 
% [1] D.L.D Caspar and A. Klug Physical Principles in the Construction of
% Regular Viruses. Cold Springs Harb. Symp. Quant. Biol. 27, 1962

%% Generate the points on the north polar triangles. For each polar triangle,
%rotate the points around one of the non polar vertices twice to generate
%each of the ten equatorial triangles

X = [];

RotNP = [cos(2*pi/5),-sin(2*pi/5),0;
         sin(2*pi/5),cos(2*pi/5),0;
         0          ,0          ,1];

Roty = [1/sqrt(5),0, -2/sqrt(5);
        0        ,1,0          ;
        2/sqrt(5),0,1/sqrt(5)  ];

  Y = TriRadialAreaProj(m,n)';    
for k=0:4
 
 Rotz = [cos(pi/5+k*2*pi/5),sin(pi/5+k*2*pi/5),0;
         -sin(pi/5+k*2*pi/5),cos(pi/5+k*2*pi/5),0;
         0                  ,0                 ,1];
    
 RotVert = Rotz^(-1)*Roty^(-1)*RotNP*Roty*Rotz;

X = [X;(RotNP^k*Y)'];
X = [X;(RotVert*RotNP^k*Y)'];
X = [X;(RotVert^2*RotNP^k*Y)'];

%To generate south polar triangles, rotate base triangle around x
%axis by pi and shift around south pole by pi/5+k*2pi/5.

Rotx = [1,0,0;
        0,-1,0;
        0,0,-1];

X = [X;(Rotz*Rotx*Y)'];    
end

%Remove duplicate points
X = real(X);
X = 10^(8).*arrayfun(@round,10^(8).*X);
X = unique(X,'rows');

%Triangulate the nodes
tri = delaunay(x);
tri = freeBoundary(TriRep(tri,x));

end
function [x] = TriRadialAreaProj(m,n)
if (n>m)
    k = n;
    n = m;
    m = k;
end

side = sqrt(4*pi/(5*sqrt(3)));
circum = side*sin(2*pi/5); 

%Change of base matrix from (m,n) coordinates to Euclidean
A = [1, 1/2;
     0, sqrt(3)/2];

 
 a = A*[0;0];
 b = A*[m;n];
 c = A*[-n;m+n];
 

%Generate lattice on plane

e1 = [1,0];
e2 = [1/2,sqrt(3)/2];

L = ones(m+n+1,1)*e1;
L=[-n:m;-n:m]'.*L;

for i = 0:m+n
         L1 = L(1:m+n+1,:) + i*ones(m+n+1,1)*e2;
         L = [L;L1];
end

 [sizeL,~] = size(L);

 Mbc = (c(2,1)-b(2,1))/(c(1,1)-b(1,1));
 
 %Change of base matrix sending [m;n] to [1;0] and [-n;m+n] to [0;1]
 MN = 1/(m^2+m*n+n^2)*[m+n,n;
      -n,m];
  

T = circum*[2/sqrt(5)*cos(-pi/5),2/sqrt(5)*cos(pi/5);
     2/sqrt(5)*sin(-pi/5),2/sqrt(5)*sin(pi/5);
      1/sqrt(5)-1,1/sqrt(5)-1];

%Rot represents rotation by 2pi/3 around the center of the triangle
ctao = (5+2*sqrt(5))/sqrt(((5+sqrt(5))^2+(5+2*sqrt(5))^2));

Roty2 = [ctao,0,-sqrt(1-ctao^2);
        0,       1, 0       ;
        sqrt(1-ctao^2),0,ctao ];

Rotz2 = [-1/2,-sqrt(3)/2,0;
        sqrt(3)/2,-1/2,0 ;
        0        ,0   ,1];

Rot = Roty2^(-1)*Rotz2*Roty2;    
RotInv = Roty2^(-1)*Rotz2^(-1)*Roty2;
    
p=1;

%For each point in the lattice, test if it is in the triangle. If so,
%rotate it and project it.
for j=1:sizeL
    %Determine which points lie in the triangle third defined by two
    %vertices and the face center. Due to precision errors in Matlab, round
    %the Boolean statements to a reasonable tolerance.
    rL = round(1e10*L(j,2));
    rMbc = round(1e10*(Mbc*(L(j,1)-b(1,1))+b(2,1)));
    rA = round(1e10*norm(L(j,:)-a'));
    rB = round(1e10*norm(L(j,:)-b'));
    rC = round(1e10*norm(L(j,:)-c'));
    if(rL <= rMbc && rA>=rC && rA>=rB)
        XIcos = [(T*MN*A^(-1)*[L(j,1);L(j,2)]+[0;0;circum])'];
        xy = XIcos(1,:)/norm(XIcos(1,:));
        x(p,:) = xy;
        x(p+1,:) = Rot*xy';
        x(p+2,:) = RotInv*xy';
        p=p+3;
        
    end
    end
end
    
