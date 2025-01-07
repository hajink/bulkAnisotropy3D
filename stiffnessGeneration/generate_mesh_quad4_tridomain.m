function [xn,ien] = generate_mesh_quad4_tridomain(nsd,nen,nel,nnp,nnx,nny,lx,ly)

% coordinates
xn = zeros(nsd,nnp);

nelx = nnx-1;
nely = nny-1;

for i = 1:nnx
    xn(1,i) = (i-1)*(lx/nelx); 
end

dy = ly/nely;
dx = dy/atan(ly/lx);

dy 

dx
tmp = lx-2*dx
hx = lx/nelx
nelx_new = 2*round(tmp/(2*hx))
de = (nelx-nelx_new)/2

nnp = nn

loc = (de:(nnx-de))+nnp;
xn(1,i) = (i-1)*(lx/nelx); 
xn(2,loc) = (i-1)*(ly/nely); 
return

for i = 2:nny
    dx = 1
    loc = (1:nnx)+(i-1)*nnx;
    xn(1,loc) = xn(1,1:nnx);
    xn(2,loc) = (i-1)*(ly/nely); 
end

% plot(xn(1,:),xn(2,:),'*')
% axis equal

% connectivity
ien=zeros(nen,nel);

for i = 1:nelx
    ien(:,i) = [0  1  1+nnx  nnx]'+i; 
end
for i = 2:nely
    loc = (1:nelx)+(i-1)*nelx;
    ien(:,loc) = ien(:,1:nelx)+(i-1)*nnx; 
end
return
