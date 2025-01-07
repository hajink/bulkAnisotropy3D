% =======================================================================

function [Ke0] = initialize_FE_simplified(var, f,idata,idb,...
    ien,iheat,iplane,isolid,l,matno,ndf,ngp,nel,nen,Nint,nnp,nsd,nstre,...
    point,props,props_pl,weight,xn)

% Initializes the FE solver
% ------------------------------------------------------------------
% f         = applied force
% idata     = problem number
% idb       = flag for prescribed boundary conditions
% ien       = element connectivity
% iheat     = identifier for conduction problems
% iplane    = 1 - plane strain, 2 - plane stress
% isolid    = identifier for solids only modeling
% l         = location array for complaint design problems
% matno     = number of solid materials
% ndf       = number of degrees of freedom
% ngp       = number of gauss points per element
% nel       = number of elements
% nen       = number of nodes per element
% nnp       = number of nodes
% nsd       = number of spacial dimensions
% nstre     = number of stresses in each integration point
% props     = elastic properties
% props_pl  = plastic properties
% xn        = nodal coordinates
%
% D         = elastic constitutive matrix
% Eelem     = Young's modulus of elements
% elem      = map of elements that each node is talking to (solids only)
% F         = global load vector
% Helem     = hardening moduli for the elements
% L         = location vector for complaint mechanism problems
% LM        = local to global mapping
% id
% Imx       = identity matrix
% mu_mat    = shear moduli of the elements
% nee       = number of element equations
% nelem     = number of elements that each node is talking to (solids only)
% neq       = number of equations
% Ke0       = solid element stiffness matrices
% P         = projection matrix
% snuelem   = Poisson's ratio for elements
% sum_nodes = how many elements each node is talking to (solids only)
% syelem    = yield stress of elements
% tgp       = total number of gauss points
% thick     = element thicknesses
%
% ------------------------------------------------------------------

% =======================================================================

% -------------------------------------------------------------------------
% NUMBER EQUATIONS
% -------------------------------------------------------------------------
% number the equations
[id,neq]=number_eq(idb,nnp,ndf);

% Organize equation number information to element level
nee_tmp = nsd*nen;
nee     = ndf*nen;   % Number of element equations
tgp     = nel*ngp;   % total number of gauss points
loc     = 1:nee;
LM  = sparse(nee_tmp,nel);
for ielem=1:nel
    [LM(loc,ielem)] = get_local_id(id,ien(:,ielem),nen,ndf);
end

% -------------------------------------------------------------------------
% ASSEMBLE GLOBAL LOAD VECTOR
% -------------------------------------------------------------------------
% compute the constitant element load vector for each element
if (neq > 0)
    [F] = loadps(f,id,neq,nnp,ndf);

    L = 0;
    % for complaint mechanism problems set up the location vector
    if idata >= 21 && idata <= 25
        [L] = loadps(l,id,neq,nnp,ndf);
    end
end


% -------------------------------------------------------------------------
% COMPUTE THE ELEMENT STIFFNESS MATRICES
% -------------------------------------------------------------------------
Ke0      = zeros(nee_tmp,nee_tmp);
thick    = sparse(nel,1);
Eelem    = sparse(nel,1);
snuelem  = sparse(nel,1);
syelem   = sparse(nel,1);
Helem    = sparse(nel,1);
mu_mat   = sparse(nel,1);
D        = zeros(nstre,nstre);
Imx_nsd  = eye(nsd);
zero_nsd = zeros(nsd);
coeff_q  = sparse(tgp,1);
kgp = 0;

% solid stiffness
for ielem = 1

    % unpack the material properties for the element
    [E,snu,t,sy0,H,lprop] = unpack_mat(ielem,props,props_pl,matno);

    % save the parameters
    thick   = t*ones(nel,1);
    Eelem   = E*ones(nel,1);
    syelem  = sy0*ones(nel,1);
    snuelem = snu*ones(nel,1);
    Helem   = H*ones(nel,1);

    % compute the shear modulus
    mu_mat    = E/(2*(1+snu))*ones(nel,1);

    % get number of stress components
    nstr = nsd*(nsd+1)/2;

        % construct the elastic constitutive matrix
        D = D_3d(var,snu,iplane,nstr);

        Ke0(loc,loc) = Ke_bric8(var,iplane,snu,nee,nen,nsd,...
            ndf,ien(:,ielem),xn,Imx_nsd,zero_nsd);
    end
% end
coeff_q = repmat(coeff_q,1,nee);


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% HELPER FUNCTION DEFINITIONS BELOW:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [id,neq]=number_eq(idb,nnp,ndf)
    
    %------------------------------------------------------------------------
    %  Assigns equation numbers to unknown displacement degrees of freedom
    %     for global stiffness and force assembly
    %
    %  ASSUMES REDUCED STORAGE APPROACH!
    %  Need more output for Complete Storage Approach
    %
    %  Variable Descriptions:
    %  Return:
    %     id(i,n) = equation number corresponding to dof i of node n
    %     neq = total number of equations
    %  Given:
    %     idb(i,n) = boundary condition flag for dof i of node n
    %     nnp  = number of nodes
    %     ndf  = number of degrees of freedom per node
    %------------------------------------------------------------------------
    
    id=zeros(ndf,nnp);  % initialize id
    neq=0;               % number of equations
    
    for n = 1:nnp
    
        % check if boundary conditions are applied
        tmp = find(idb(1:ndf,n)==0);
    
        for i = 1:length(tmp)
    
            % udate # of equations
            neq = neq+1;
    
            % give equation number
            id(tmp(i),n) = neq;
        end
    end
    
    
    % for n = 1:nnp
%     for i = 1:ndf
%         if idb(i,n) == 0
%             % udate # of equations
%             neq = neq + 1;
%
%             % if no prescribed displacement at dof i of node n
%             %   => give an equation # different from 0
%             id(i,n) = neq;
%
%         end
%     end
% end


%********************************************************************%
function [LM]=get_local_id(id,ien,nen,ndf)
    %------------------------------------------------------------------------
    %  Localizes the global equations numbers for a single element
    %     for easier global stiffness and force assembly
    %     in main: LM(n,e)=global eqn number where n=local eqn number, e=element
    %
    %  Variable Descriptions:
    %  Return:
    %     LM(n) = global eqn number corresponding to local eqn number n
    %  Given:
    %     id(i,n) = global equation number dof i of global node n
    %     ien(i,e) = global node number for local node i of elment e
    %     nen  = number of nodes per element
    %     ndf  = number of degrees of freedom per node
    %------------------------------------------------------------------------
    
    nee = ndf*nen;              % # of element equations
    LM  = zeros(nee,1);         % initialize LM
    
    k = 0;                      % counting local equation #s.
    
    % Global equation numbers for the element
    for i = 1:nen                       % Loop over # of nodes per element
        for j = 1:ndf                   % Loop over # of dof's per node
            k = k + 1;
            LM(k) = [ id(j,ien(i)) ];
        end
end

%****************************************************
function[F] = loadps(f,id,neq,nnp,ndf)
    % function that defines the global  nodal loads
    
    % initiate global load vecctor
    F = sparse(neq,1);
    
    % read nodal loads
    F = add_loads_to_force(F,f,id,nnp,ndf);
    
    % % Compute forces from ds ~= 0 and insert into F
    % Fse = zeros(nee,nel);  dse = zeros(nee,nel);  % We will store all nel element vectors
    % for i=1:nel
    %     dse(:,i) = get_de_from_dcomp(g,ien(:,i),nen,ndf);  % get dse for current element
    %     Fse(:,i) = - Ke(:,:,i) * dse(:,i);                 % compute element force
    %     F = addforce(F,Fse(:,i),LM(:,i),nee);     % assemble elem force into global force vector
    % end
    % clear Fse dse;

return

%********************************************************************************
function [F]=add_loads_to_force(F,fn,id,nnp,ndf)
    
    %------------------------------------------------------------------------
    %  Adds the applied load nodal vector into global force vector using eqn numbers
    %
    %  Variable Descriptions:
    %  Return:
    %     F  = updated global force vector
    %  Given:
    %     F  = global force vector
    %     fn = applied nodal loads (ndf,nnp)
    %     id(i,n) = displacement equation number for dof i, node n 
    %     nnp = number of nodal points
    %     ndf = number degrees of freedom per node
    %------------------------------------------------------------------------
    
    % Loop over nodes and degrees of freedom
    for n=1:nnp
        for i=1:ndf
            M = id(i,n);      % Get Global equation number 
            if (M > 0)        % Check Global Equation Number
                F(M) = F(M) + fn(i,n);
            end
        end
    end
return

%****************************************************
function [E,snu,thick,sy0,H,lprop] = unpack_mat(ielem,props,plastic_prop,matno)
    
    lprop = matno(ielem);
    E     = props(lprop,1);
    snu   = props(lprop,2);
    thick = props(lprop,3);
    
    if size(plastic_prop,1)>0
        sy0 = plastic_prop(:,1,lprop);
        H   = plastic_prop(:,2,lprop);
    else
        sy0 = 0;
        H   = 0;
    end
return

%********************************************************************************
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

function [kee] = Ke_bric8(var,iplane,snu,nee,nen,nsd,ndf,ien,xn,Imx_nsd,...
    zero_nsd)

    % =======================================================================
    
    % Define location and weights of integration points
    [Nint,point,weight] = gausspoints_quad4(ndf,nsd);
    
    % get constitutive matrix
    nstr = nsd*(nsd+1)/2;
    D = D_3d(var,snu,iplane,nstr);
    
    % Compute element matrix
    [kee] = ke_elem_bric8(Nint,point,weight,ien,xn,ndf,nen,nee,D,Imx_nsd,...
        zero_nsd);

return

% =======================================================================

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
return 

%*****************************************************************************