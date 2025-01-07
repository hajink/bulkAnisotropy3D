% Code for shape functions N for Bric8 elements

% =======================================================================

function [sN,dNdx,dNdy,dNdz,detJ,invJ] = shape_bric8(q,r,s,xn,ien,nen,...
    Imx_nsd,zero_nsd)

% =======================================================================
%  Computes the shape functions for quad4 elements
% =======================================================================

% Initialize
dNdx  = zeros(nen,1);
dNdy  = dNdx;
dNdz  = dNdx;

% compute shape functions at(r,s)
[sN,dNdq,dNdr,dNds] = shapefunct_bric8(q,r,s);

% compute Jacobian, its deterninant and inverse
[~,detJ,invJ] = jacobian_3d(dNdq,dNdr,dNds,xn(:,ien),nen,Imx_nsd,zero_nsd);

% compute the derivatives of the shape functions
for i = 1:nen
    dNdx(i) = invJ(1,1)*dNdq(i)+invJ(1,2)*dNdr(i)+invJ(1,3)*dNds(i);
    dNdy(i) = invJ(2,1)*dNdq(i)+invJ(2,2)*dNdr(i)+invJ(2,3)*dNds(i);
    dNdz(i) = invJ(3,1)*dNdq(i)+invJ(3,2)*dNdr(i)+invJ(3,3)*dNds(i);
end

return
      
% =======================================================================







