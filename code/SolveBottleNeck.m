% Returns the voltage and electric field components for a bottle neck
% simulation
function [Vmap, Ex, Ey] = SolveBottleNeck(Xmax,Ymax,nx,ny,sigma1,sigma2,Vapp,xlimleft, xlimright, ylimtop, ylimbot )

sMap = ones(nx,ny)*sigma1;

for i = 1:nx
    x = i/nx*Xmax;
    for j = 1:ny
        y = j/ny*Ymax;
        if (x >= xlimleft && x <= xlimright && (y >= ylimtop || y <= ylimbot))
            sMap(i,j) = sigma2;
        end
    end
end

G = zeros(nx*ny);
B = zeros(1,nx*ny);

for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        if i == 1
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = Vapp;
        elseif i == nx
            G(n,:) = 0;
            G(n,n) = 1;
        elseif j == 1
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (sMap(i,j) + sMap(i-1,j))/1.0;
            rxp = (sMap(i,j) + sMap(i+1,j))/1.0;
            ryp = (sMap(i,j) + sMap(i,j+1))/2.0;

            G(n,n) = -(rxm+rxp+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
        elseif j ==  ny
            
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;

            rxm = (sMap(i,j) + sMap(i-1,j))/1.0;
            rxp = (sMap(i,j) + sMap(i+1,j))/1.0;
            rym = (sMap(i,j) + sMap(i,j-1))/2.0;

            G(n,n) = -(rxm+rxp+rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            
            rxm = (sMap(i,j) + sMap(i-1,j))/2.0;
            rxp = (sMap(i,j) + sMap(i+1,j))/2.0;
            rym = (sMap(i,j) + sMap(i,j-1))/2.0;
            ryp = (sMap(i,j) + sMap(i,j+1))/2.0;
            
            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
        end        
    end
end
G = sparse(G);
V = G\B';

Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        Vmap(i,j) = V(n);
    end
end

dx = Xmax/nx;
dy = Ymax/ny;
Ex = zeros(nx,ny);
Ey = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i,j) = (Vmap(i+1,j) - Vmap(i,j));
        elseif i == nx
            Ex(i,j) = (Vmap(i,j) - Vmap(i-1,j));
        else
            Ex(i,j) = (Vmap(i+1,j) - Vmap(i-1,j))*0.5;
        end
        if j == 1
            Ey(i,j) = (Vmap(i,j+1) - Vmap(i,j));
        elseif j == ny
            Ey(i,j) = (Vmap(i,j) - Vmap(i,j-1));
        else
            Ey(i,j) = (Vmap(i,j+1) - Vmap(i,j-1))*0.5;
        end
    end
end

Ex = sMap.*(-Ex/dx);
Ey = sMap.*(-Ey/dy);

end