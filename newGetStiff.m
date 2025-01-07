function[Ke0] = newGetStiff(var, nelx, nely, nelz)

%Cut from topopt file in JVC code. Just need Ke0

idata = 5;
% unpack problem data
[f,idb,ien,iheat,l,matno,ndf,nel,nen,nnp,nsd,props,props_pl,xn] = benchmark_problems(idata,var); %currently defines lx = nelx;

props(1,1:2) = [var.E1 var.nu]; %replace properties with desired problem data

x = repmat(var.volfrac,[nely,nelx,nelz]);
isolid = 0;
iplane = var.isotropy;
    % get gauss point information
    [Nint,point,weight,nstre(1),nstr1,ngp] = gausspoints_quad4(ndf,nsd);

        [Ke0] = initialize_FE_simplified(var,f{1},idata,...
            idb{1},ien{1},iheat(1),iplane,isolid,l,matno,ndf(1),ngp,nel,...
            nen,Nint,nnp,nsd,nstre(1),point,props(1,:),props_pl,weight,xn);

return