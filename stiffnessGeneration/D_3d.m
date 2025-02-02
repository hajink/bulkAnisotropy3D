function[D] = D_3d(var,snu,iplane,nstr)

% initialize and define D
D = zeros(nstr,nstr);


if (iplane == 0)         % isotropic    
    coeff      = var.E1/((1+snu)*(1-2*snu));
    D(1:3,1:3) = snu;
    D(1,1)     = 1-snu;
    D(2,2)     = D(1,1);
    D(3,3)     = D(1,1);
    D(4,4)     = 0.5*(1-2*snu);
    D(5,5)     = D(4,4);
    D(6,6)     = D(4,4);
    D          = coeff*D;
end
if (iplane == 1)        % anisotropic
    %make new D with multiple E values (1, 2, or 3)
    E1 = var.E1;  %x
    E2 = var.E2;   %y
    E3 = var.E3;   %z
    
    G1 = E1/(2*(1+snu));
    G2 = E2/(2*(1+snu));
    G3 = E3/(2*(1+snu));
    
    coeff = (1-snu*snu - snu*snu - snu*snu - 2*snu*snu*snu)/(E1*E2*E3);
    D(1,1) = (1-snu*snu)/(E2*E3*coeff);
    D(1,2:3) = (snu+snu*snu)/(E2*E3*coeff);
    D(2,1) = (snu+snu*snu)/(E3*E1*coeff);
    D(2,2) = (1-snu*snu)/(E3*E1*coeff);
    D(2,3) = D(2,1);
    D(3,1:2) = (snu+snu*snu)/(E1*E2*coeff);
    D(3,3) = (1-snu*snu)/(E1*E2*coeff);
    D(4,4)= 2*G1;
    D(5,5) = 2*G2;
    D(6,6) = 2*G3;
end

return