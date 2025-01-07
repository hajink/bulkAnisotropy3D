% =======================================================================

function [Jaco,detJ,invJ] = jacobian_3d(dNdq,dNdr,dNds,coord,nen,Imx,Jaco)

% =======================================================================
%  Computes the Jacobian, its determinate and inverse
% =======================================================================

% Add to Jacobian
for i = 1:nen
    Jaco(1,1) = Jaco(1,1)+dNdq(i)*coord(1,i);
    Jaco(1,2) = Jaco(1,2)+dNdq(i)*coord(2,i);
    Jaco(1,3) = Jaco(1,3)+dNdq(i)*coord(3,i);

    Jaco(2,1) = Jaco(2,1)+dNdr(i)*coord(1,i);
    Jaco(2,2) = Jaco(2,2)+dNdr(i)*coord(2,i);
    Jaco(2,3) = Jaco(2,3)+dNdr(i)*coord(3,i);

    Jaco(3,1) = Jaco(3,1)+dNds(i)*coord(1,i);
    Jaco(3,2) = Jaco(3,2)+dNds(i)*coord(2,i);
    Jaco(3,3) = Jaco(3,3)+dNds(i)*coord(3,i);
    
    
    % Jaco(1,1) = Jaco(1,1)+dNdq(i)*coord(1,i);
    % Jaco(1,2) = Jaco(1,2)+dNdr(i)*coord(1,i);
    % Jaco(1,3) = Jaco(1,3)+dNds(i)*coord(1,i);
    % 
    % Jaco(2,1) = Jaco(2,1)+dNdq(i)*coord(2,i);
    % Jaco(2,2) = Jaco(2,2)+dNdr(i)*coord(2,i);
    % Jaco(2,3) = Jaco(2,3)+dNds(i)*coord(2,i);
    % 
    % Jaco(3,1) = Jaco(3,1)+dNdq(i)*coord(3,i);
    % Jaco(3,2) = Jaco(3,2)+dNdr(i)*coord(3,i);
    % Jaco(3,3) = Jaco(3,3)+dNds(i)*coord(3,i);
end

% Find its determinant
detJ = det(Jaco);

% Find its inverse
invJ = Imx/Jaco;

return

% =======================================================================