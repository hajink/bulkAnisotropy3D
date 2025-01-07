% =======================================================================

function[f,idb,ien,iheat,l,matno,ndf,nel,nen,nnp,nsd,props,props_pl,xn] = benchmark_problems(idata,var)

% Sets up benchmark problems
% ------------------------------------------------------------------
% idata     = problem number
% nelx      = number of elements in x
% nely      = number of elements in y
% volfrac   = specified volume fraction
%
% dr        = minimum length scale deviation for robust problems
% error_tol = error tolerance for symmetry constraints of homogenization
%             problems
% f         = applied force
% fac_obj   = factor for objective function
% g         = magnitude of prescribed boundary conditions
% idb       = flag for prescribed boundary conditions
% idispbc   = flag for displacement controlled problem
% ien       = element connectivity
% ien_dv    = element connectivity for design variable plot
% iobj      = objective marker
% iiobj     = objective for elastic material design problem
% iheat     = conductive problem identifier
% ihom      = homogenization identifier
% imultmat  = multimaterial identifier
% iplane    = 1 - plane strain, 2 - plane stress
% ipre_rho  = identifier for prescribed design variables
% irobust   = robust formulation identifier
% isolver   = type of solver
% isym      = symmetry constraints for homogenization problems
% ivol      = identifier for volume constraint types for multimaterial
%             problems
% l         = locator of displacements to be minmized for inverter problems
% lx        = length of design domain in x
% ly        = length of design domain in y
% matno     = number of solid materials
% ncase     = number of test strain cases for homogenization problems
% ncells    = number of cells in periodic samples
% ndf       = number of degrees of freedom
% ndv       = number of design variables
% ndv_x     = number of design variables in x
% ndv_y     = number of design variables in y
% ndv_z     = number of design variables in z
% nel       = number of elements
% nen       = number of nodes per element
% nnlat     = lattice for element print
% nnlat_dv  = lattice for design variable print
% nnp       = number of nodes
% nnp_dv    = number of nodes for design variable plot
% nnx       = number of nodes in x
% nny       = number of nodes in y
% nobj      = number of objectives
% nsd       = number of spacial dimensions
% props     = elastic properties
% props_pl  = plastic properties
% spring    = spring locations and constants for inverter problems
% V         = volume constraint values
% v_elem    = element volume
% v_tol     = tolerance on volume constraint |V_actual - V_goal| < v_tol
% xn        = nodal coordinates
% xn_dv     = nodal coordinates for design variable plot
%
% ------------------------------------------------------------------

% =======================================================================
volfrac = var.volfrac;
nelx = var.nelx;
[nely, nelz] = var.getDims();
rmin = var.rmin;
nsd =  3;          % number of spacial dimensions
ndf =  nsd;        % number of degrees of freedom per node
nen =  4*(nsd-1);  % number of element nodes
ngp =  4*(nsd-1);  % number of gauss points per element

nelz = (nsd-2)*nelz+1-mod(nsd,2);
nel  = nelx*nely*nelz; % # of elements
nnx  = nelx+1;  % number of nodes in x direction
nny  = nely+1;  % number of nodes in y direction
nnz  = nelz+1;  % number of nodes in z direction
nnp  = nnx*nny; % number of nodal pts

ndv_x = nnx;   % number of design variables in x
ndv_y = nny;   % number of design variables in y
ndv_z = 0;     % number of design variables in z

% update if 3d
if nsd == 3
    ndv_z = nnz;   
    nnp = nnp*nnz;
end

% set number of design variables = number of nodes
ndv   = nnp;

% number of load cases
ncase = 1;

isym_geo  = 0; % 0 - no geometric symmetry, 1 - geometric symmetry enforced

% location vector for complaint design problems and stiffness constants
l      = 0;

% conductive problem identifier
iheat = 0;

% choose periodic sample size
nperx = 1;
npery = nperx;

% projection identifier
% min length scale applied on 0 - the solid phase only, 1 - the void phase
% only, 3 - both solid and void phases, 6 - discrete objects
MLS = 0;


    [f,g,h,idb,ien,lx,ly,lz,matno,ndv,nel,nen,props,...
        props_pl,V,v_elem,xn] = benchmark_iobj5(idata,MLS,ncase,ndf,ndv,nen,...
        nel,nelx,nely,nelz,nnp,nnx,nny,nnz,nsd,rmin,volfrac);
    nel_t = nel;
    

if isym_geo > 0
    if nsd == 2
        ndv_x = floor((ndv_x-1)/2+1);
        ndv_y = floor((ndv_y-1)/2+1);
        ndv = ndv_x*ndv_y;
    end
end

% generate mesh for design variable plots
if ndv == nnp && nen(1) > 2
    nnx_dv = nnx+1;
    nny_dv = nny+1;
    nnz_dv = nnz+1;
    if nsd == 2
        nnp_dv = nnx_dv*nny_dv;
        [xn_dv,ien_dv] = generate_mesh_quad4(nsd,nen(1),nnp,nnp_dv,...
            nnx_dv,nny_dv,lx+h,ly+h);
    else
        nnp_dv = nnx_dv*nny_dv*nnz_dv;
        [xn_dv,ien_dv] = generate_mesh_bric8(nsd,nen(1),nnp,nnp_dv,...
            nnx_dv,nny_dv,nnz_dv,lx+h,ly+h,lz+h);
    end
else
    xn_dv  = xn;
    ien_dv = ien;
    nnp_dv = nnp;
end

% prepare lattices for vtk files
nnlat    = [nnx  nny  nnz]';
return


function[f,g,h,idb,ien,lx,ly,lz,matno,ndv,nel,nen,props,props_pl,V,...
    v_elem,xn] = benchmark_iobj5(idata,MLS,ncase,ndf,ndv,nen,nel,nelx,nely,nelz,nnp,...
    nnx,nny,nnz,nsd,rmin,volfrac)

% half mbb beam

% in mm
lx = nelx; %HAJIN OVERRIDE
ly = nely;
lz = nelz;

% if nelx/nely ~= lx/ly
%     disp('error in mesh dimensions')
%     return
% end

if idata == 5
    
    % generate mesh
    if nsd == 2
        [xn,ien_tmp] = generate_mesh_quad4(nsd,nen,nel,nnp,nnx,nny,lx,ly);
        loc_load = nely*nnx+1;
    else
        %HAJIN changed mesh generation (from generate_mesh_bric8 to generate_mesh_top3d):
        [xn,ien_tmp] = generate_mesh_top3d(nsd,nen,nel,nnp,nnx,nny,nnz,...
            lx,ly,lz);
        % [xn,ien_tmp] = generate_mesh_bric8(nsd,nen,nel,nnp,nnx,nny,nnz,...
        %     lx,ly,lz);
        loc_load = nely*nnx+1:nnx*nny:nnp;
    end
    ien{1} = ien_tmp;  
else
    
    % truss mesh
    nen = 2; % number of nodes per element
    
    itruss = 7;     % 0 - local mesh for 2D
    % 1 - global connetivity,
    % 2 - local connectivity with no diagonals,
    % 3 - local connectivity with internal diagonals,
    % 4 - local connectivity with face diagonals,
    % 5 - local connectivity with all diagonals,
    % 6 - local connectivity with only diagonals,
    % 7 - local connectivity within adjecent cells   
    
    % gernerate mesh and connectivity
    [xn,ien,nel,nnp] = groundTruss(nnx,nny,1,nen,nsd,lx,ly,1,itruss);
    
    ndv      = nel;
    DeltaE   = ones(nel,1);
end

% define materials

nomats = 1;             % number of materials

% material 1 elatsic properties THESE GET OVERRIDDEN 
E      = 1;             % Young's modulus
snu    = 0.30;          % Poisson's ratio
if idata == 105
    snu = 10; % truss bar area
end

if MLS == 3
    E = E*100;
end

% geometric properties
t      = 1;
h      = lx/nelx;
v_elem = (lx/nelx)*(ly/nely)*(lz/nelz)*ones(nel,1);
V      = volfrac*lx*ly*lz;

% assign prooerties to elements
for i = 1:nomats
    props(i,:) = [E(i)  snu(i)  t(i)  0];
end
props_pl = [];
matno = ones(nel,1);

% define boundary conditions and applied force (symmetry and fixed at end)
idb1                       = zeros(ndf(1),nnp,ncase);
g1                         = zeros(ndf(1),nnp,ncase);
f1                         = zeros(ndf(1),nnp,ncase);
idb1(1,1:nnx:nnp)          = 1;
idb1(2,(nnx-round(rmin/h)+1):nnx*nny:nnp)    = 1;
f1(2,loc_load)             = -1;
if nsd == 3
    idb1(3,(nnp-nnx*nny):nnp) = 1;
    ipre_rho = 65;
end

% set prescribed rho if relevant
ipre_rho = 65;

% save idb, g, f as cells
idb = cell(ncase,1);
g   = cell(ncase,1);
f   = cell(ncase,1);

idb{1} = idb1;
f{1} = f1;
g{1} = g1;
return
