% =======================================================================

function [kee] = ke_elem_bric8(Nint,point,weight,ien,xn,ndf,nen,nee,D,...
    Imx_nsd,zero_nsd)

% =======================================================================


% Initialize 
kee   = zeros(nee,nee);
coeff = 0;

for i = 1:Nint(1)
    qi = point(i,1);
    wi = weight(i,1);
    
    for j = 1:Nint(2)
        rj = point(j,2);
        wj = weight(j,2);
        
        for k = 1:Nint(3)
            sk = point(k,3);
            wk = weight(k,3);
            
            [~,dNdx,dNdy,dNdz,detJ]  = shape_bric8(qi,rj,sk,xn,ien,nen,...
                Imx_nsd,zero_nsd);
            
            B = B_3d_elastic(dNdx,dNdy,dNdz,nen,ndf);
                  
            coeff = wi*wj*wk*detJ;
            
            % Compute kee = kee + transpose(B)*D*B*coeff
            kee = kee+(B')*D*B*coeff;
        end
    end
end

end
      
% ======================================================================= 