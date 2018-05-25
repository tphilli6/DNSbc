
Nx=1000;
Ny=65;
Nz=65;

x0=0;
dt=0.003;
x=0:dt:Nx*dt;

y=cos( (Ny-1:-1:0)*pi/(Ny-1) );

z0=0;
zL=4/6*pi;
z=linspace(z0,zL,Nz);

[Y,Z]=meshgrid(y,z); Y=Y'; Z=Z';
[YY, ZZ, XX] = meshgrid(y, z, x);

data=importdata('fort.77');

D=reshape(data,[Ny,Nz,Nx]);
% clear data

% for i=1:Nx
%    surf(Z, Y, D(:,:,i) )
%    view([0,0,90])
%    shading interp
%    
%    pause(0.01)
% end

Lzcorrsum=zeros(Ny,Nz-1);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz-1
          Lzcorr(k)= sum(D(j,k:end,i).*D(j,1:end-k+1,i) );
        end
        Lzcorr=Lzcorr/Lzcorr(1);
        Lzcorrsum(j,:) = Lzcorrsum(j,:) + Lzcorr;

    end
end
Lzcorrsum = Lzcorrsum/Nx;

for j=1:Ny
Lz(j) = interp1(Lzcorrsum(j,:), z(1:end-1), exp(-pi/4));
end