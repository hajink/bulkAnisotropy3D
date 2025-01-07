% =======================================================================

function [sN,dNdq,dNdr,dNds] = shapefunct_bric8(q,r,s)

% =======================================================================
%  Evaluates the shape functions and their derivatives at point (q,s,r)
% =======================================================================

% Initialize
sN   = zeros(8,1);
dNdq = sN;
dNdr = sN;
dNds = sN;

% Shape functions
sN(1) = (1-q)*(1-r)*(1-s);
sN(2) = (1+q)*(1-r)*(1-s);
sN(3) = (1+q)*(1+r)*(1-s);
sN(4) = (1-q)*(1+r)*(1-s);
sN(5) = (1-q)*(1-r)*(1+s);
sN(6) = (1+q)*(1-r)*(1+s);
sN(7) = (1+q)*(1+r)*(1+s);
sN(8) = (1-q)*(1+r)*(1+s);
sN    = 0.125*sN;

% Derivatives
dNdq(1) = -(1-r)*(1-s);
dNdq(2) =  (1-r)*(1-s);
dNdq(3) =  (1+r)*(1-s);
dNdq(4) = -(1+r)*(1-s);
dNdq(5) = -(1-r)*(1+s);
dNdq(6) =  (1-r)*(1+s);
dNdq(7) =  (1+r)*(1+s);
dNdq(8) = -(1+r)*(1+s);
dNdq    =  0.125*dNdq;

dNdr(1) = -(1-q)*(1-s);
dNdr(2) = -(1+q)*(1-s);
dNdr(3) =  (1+q)*(1-s);
dNdr(4) =  (1-q)*(1-s);
dNdr(5) = -(1-q)*(1+s);
dNdr(6) = -(1+q)*(1+s);
dNdr(7) =  (1+q)*(1+s);
dNdr(8) =  (1-q)*(1+s);
dNdr    =  0.125*dNdr;

dNds(1) = -(1-q)*(1-r);
dNds(2) = -(1+q)*(1-r);
dNds(3) = -(1+q)*(1+r);
dNds(4) = -(1-q)*(1+r);
dNds(5) =  (1-q)*(1-r);
dNds(6) =  (1+q)*(1-r);
dNds(7) =  (1+q)*(1+r);
dNds(8) =  (1-q)*(1+r);
dNds    =  0.125*dNds;

return

% =======================================================================