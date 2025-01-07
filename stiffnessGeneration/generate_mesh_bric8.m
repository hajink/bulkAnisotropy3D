function [xn,ien] = generate_mesh_bric8(nsd,nen,nel,nnp,nnx,nny,nnz,lx,ly,lz)

% coordinates
xn = zeros(nsd,nnp);

nelx = nnx-1;
nely = nny-1;
nelz = nnz-1;

for i = 1:nnx
    xn(1,i) = (i-1)*(lx/nelx); 
end
for i = 2:nny
    loc = (1:nnx)+(i-1)*nnx;
    xn(1,loc) = xn(1,1:nnx);
    xn(2,loc) = (i-1)*(ly/nely); 
end
tmp = nnx*nny;
for i = 2:nnz
    loc = (1:tmp)+(i-1)*tmp;
    xn(1:2,loc) = xn(1:2,1:tmp); 
    xn(3,loc)   = (i-1)*(lz/nelz); 
end

% connectivity
ien = zeros(nen,nel);
for i = 1:nelx
    ien(1:4,i) = [0  1  1+nnx  nnx]'+i;
    ien(5:8,i) = ien(1:4,i)+tmp;
end
for i = 2:nely
    loc = (1:nelx)+(i-1)*nelx;
    ien(:,loc) = ien(:,1:nelx)+(i-1)*nnx; 
end
tmp2 = nelx*nely;
for i = 2:nelz
    loc = (1:tmp2)+(i-1)*tmp2;
    ien(:,loc) = ien(:,1:tmp2)+(i-1)*tmp;
end


return
