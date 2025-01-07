function [Nint,point,weight,nstre,nstr1,ngp] = gausspoints_quad4(ndf,nsd)
% function that defines integration points and weigths for quad4 elements

Nint(1) = 2; % number of gauss points in x
Nint(2) = 2; % number of gauss points in y
Nint(3) = nsd-1;

% total number of gauss points
ngp = Nint(1)*Nint(2)*Nint(3);

% get gauss points
if nsd < 3
    [point,weight]=gauss(Nint(1),Nint(2),0,nsd);
    
else
    [point,weight]=gauss(Nint(1),Nint(2),Nint(3),nsd);
end

% define number of individual stress/strain components
nstre = nsd*(nsd+1)/2;

% define number of stress components in each element (4 in 2d, nstre in 3d)
nstr1 = nstre+1-mod(nsd,2);