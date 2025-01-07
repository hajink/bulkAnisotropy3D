% AN 169 LINE 3D TOPOLOGY OPITMIZATION CODE BY LIU AND TOVAR (JUL 2013)
% modified by Hajin Kim-Tackowiak:
% -MMA
% -Heaviside filter applied
% -eta continuation (SIMP)
% -passive region
% -New definition of stiffness matrix to apply basic bulk anisotropy (see
% newGetStiff.m function)
% -parameters are defined and contained in object class userSettings that
% is the sole input. Use run_top3d() to iterate and change settings before
% each call to top3d_anisotropy()

%note: 1 element is length 1. If you want nelx != lx, must redefine
%(probably have to redefine benchmark_problems() too)

function top3d_anisotropy(var)
addpath('MMA')
addpath('stiffnessGeneration')

% USER-DEFINED LOOP PARAMETERS
maxloop = var.maxloop;    % Maximum number of iterations
tolx = var.tolx;      % Terminarion criterion
displayflag = var.displayflag;  % Display structure flag
% USER-DEFINED MATERIAL PROPERTIES
E = var.E1;           % Young's modulus of solid material
%assuming x-dir is E0, while y and z dir are E0/3
Emin = var.Emin;      % Young's modulus of void-like material
nu = var.nu;         % Poisson's ratio
rmin = var.rmin;
volfrac = var.volfrac;
isotropy = var.isotropy;

nelx = var.nelx; [nely, nelz] = var.getDims();
tFlange = var.tFlange; tWeb = var.tWeb;

%% INITIALIZE PENALIZATION FACTORS/PARAMETERS
eta = var.eta;
eta_max = var.eta_max;
beta = var.beta; %Highly sensitive to starting beta value
beta_max = var.beta_max;
objHistory = [];

%% USER-DEFINED LOAD DOFs (MBB)
[il,jl,kl] = meshgrid(0, nely, 0:nelz);                 % Coordinates
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 1;                             % DOFs
%% USER-DEFINED SUPPORT FIXED DOFs
[iif,jf,kf] = meshgrid(nelx,0,0:nelz);                  % Coordinates
fixednid = iif*(nely+1) + kf*(nelx+1)*(nely+1)+(nely+1-jf); % Node IDs
fixeddof_1 = [3*fixednid(:)-1, 3*fixednid(:)-2]; % DOFs %free to move in x

%side face rollers
[iif2,jf2,kf2] = meshgrid(0,0:nely,0:nelz);
fixednid = iif2*(nely+1) +(nely+1-jf2) + kf2*(nelx+1)*(nely+1);
fixeddof_2 = [3*fixednid(:), 3*fixednid(:)-2]; %free to move in z

fixeddof = union(fixeddof_1, fixeddof_2);
fixeddof = sort(fixeddof);
%% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,-1,ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
%--------new KE matrix setup that can do anisotropy
KE = (newGetStiff(var, nelx, nely, nelz));
%-----------
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
%% PREPARE FILTER
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);

%% INITIALIZE ITERATION
x = repmat(volfrac,[nely,nelx,nelz]);

%define void passive region
% within list passive() elements that ==1 will be voids 
%define separate one for solids
isTherePassive = var.isPassive;
midZ = round(nelz/2);
if isTherePassive == 1
    for elz = 1:nelz
        for ely = 1:nely
            if ely > round(tFlange) && ely < (nely-tFlange) && ...
                    (elz > (midZ + round(tWeb/2)) || elz < (midZ - round(tWeb/2)))
                passive(ely,elz) = 1;
            else
                passive(ely,elz) = 0;
            end
        end
    end
else
    passive = zeros(nely,nelz);
end
passive = repmat(passive,[nelx,1,1]);
x(find(passive)) = 0;

xTilde = x;
xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta); %Heaviside: 0-1 design 
loop = 0;
loopbeta = 0;
change = 1;

%% INITIALIZE MMA OPTIMIZER
m     = 1;                % The number of general constraints.
n     = nele;             % The number of design variables x_j.
xmin  = zeros(n,1);       % Column vector with the lower bounds for the variables x_j.
xmax  = ones(n,1);        % Column vector with the upper bounds for the variables x_j.
xold1 = x(:);             % xval, one iteration ago (provided that iter>1).
xold2 = x(:);             % xval, two iterations ago (provided that iter>2).
low   = ones(n,1);        % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = ones(n,1);        % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                % The constants a_0 in the term a_0*z.
a     = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 10000*ones(m,1);  % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.

% START ITERATION
while change > tolx && loop < maxloop
    loop = loop+1;
    loopbeta = loopbeta+1;
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^eta),24*24*nele,1); %removed (E0-Emin) factor
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);

    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
    c = sum(sum(sum((Emin+xPhys.^eta).*ce)));
    dc = -eta*xPhys.^(eta-1).*ce;

    dv = ones(nely,nelx,nelz);
    %% FILTERING AND MODIFICATION OF SENSITIVITIES
    dx = beta*exp(-beta*xTilde)+exp(-beta); %Heaviside
    dc(:) = H*(dc(:).*dx(:)./Hs);
    dv(:) = H*(dv(:).*dx(:)./Hs);
    %% METHOD OF MOVING ASYMPTOTES
    xval  = x(:);
    f0val = c;
    df0dx = dc(:);
    fval  = sum(xPhys(:))/(volfrac*nele) - 1;
    dfdx  = dv(:)' / (volfrac*nele);
    [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low,upp] = ...
        mmasub(m, n, loop, xval, xmin, xmax, xold1, xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d,beta);
    %% Update MMA Variables
    xnew     = reshape(xmma,nely,nelx,nelz);
    if isTherePassive == 1
        xnew(find(passive)) = 0;
    end
    xTilde(:) = (H*xnew(:))./Hs;
    xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta); %Heaviside: 0-1 design
    xold2    = xold1(:);
    xold1    = x(:);

    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    % PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,mean(xPhys(:)),change);
    objHistory = [objHistory c];
    % PLOT DENSITIES
    if displayflag && mod(loop,10)==1, clf; display_3D(xPhys); end %#ok<UNRCH>

        %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
    if eta < eta_max && (loopbeta >= 50 || change <= 0.001)
        % beta = 2*beta;
        % d_beta = 1.1;
        % beta = min(d_beta^loop + beta, beta_max);
        %no beta continuation for MBB and Cantilever problems
        loopbeta = 0;
        change = 1;
        eta = eta + 1;
        fprintf('Parameter eta increased to %g.\n',eta);
        fprintf('Parameter beta increased to %g.\n',beta);
    end
end
clf; display_3D(xPhys);
%figure(2); plot(objHistory);

filename = var.getFilename(); 
save(filename);
end

% === DISPLAY 3D TOPOLOGY (ISO-VIEW) ===
function display_3D(rho)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.5)  % User-defined display density threshold
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]);
                hold on;
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end
% =========================================================================
% === This code was written by K Liu and A Tovar, Dept. of Mechanical   ===
% === Engineering, Indiana University-Purdue University Indianapolis,   ===
% === Indiana, United States of America                                 ===
% === ----------------------------------------------------------------- ===
% === Please send your suggestions and comments to: kailiu@iupui.edu    ===
% === ----------------------------------------------------------------- ===
% === The code is intended for educational purposes, and the details    ===
% === and extensions can be found in the paper:                         ===
% === K. Liu and A. Tovar, "An efficient 3D topology optimization code  ===
% === written in Matlab", Struct Multidisc Optim, 50(6): 1175-1196, 2014, =
% === doi:10.1007/s00158-014-1107-x                                     ===
% === ----------------------------------------------------------------- ===
% === The code as well as an uncorrected version of the paper can be    ===
% === downloaded from the website: http://www.top3dapp.com/             ===
% === ----------------------------------------------------------------- ===
% === Disclaimer:                                                       ===
% === The authors reserves all rights for the program.                  ===
% === The code may be distributed and used for educational purposes.    ===
% === The authors do not guarantee that the code is free from errors, a