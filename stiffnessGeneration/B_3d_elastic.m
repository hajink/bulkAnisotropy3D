% Code for strain displacement matrix B


% =======================================================================

function [B] = B_3d_elastic(dNdx,dNdy,dNdz,nen,ndf)

% =======================================================================
%  Computes the strain-displacement matrix for an elastic 3d problem
% =======================================================================

% Initialize
B = zeros(6,nen*ndf);

% % Compute B ORIGINAL
% for i = 1:nen
%     loc = ((i-1)*ndf+1):(ndf*i);
%     B(:,loc) = [ dNdx(i)       0        0
%                       0   dNdy(i)       0
%                       0        0   dNdz(i)
%                       0   dNdz(i)  dNdy(i)
%                  dNdz(i)       0   dNdx(i)
%                  dNdy(i)  dNdx(i)       0 ];
% end 

% Compute B REARRANGED TO MATCH TOP3D PAPER NOTATION
for i = 1:nen
    loc = ((i-1)*ndf+1):(ndf*i);
    B(:,loc) = [ dNdx(i)       0        0
                      0   dNdy(i)       0
                      0        0   dNdz(i)
                 dNdy(i)  dNdx(i)       0 
                      0   dNdz(i)  dNdy(i)
                 dNdz(i)       0   dNdx(i)
                 ];
end 

end
% =======================================================================
