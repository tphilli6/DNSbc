%% DNS BC test program
clc
clear all

niter=8000;

% r=[1.0, 1.01, 1.025, 1.05, 1.1]
% pwr=[1.9, 1.86, 1.74, 1.465, 1.23]
%
%mz=121, niter=8000;
% r=[1.0, 1.01]
% pwr = [1.95, 1.60]
pwr=1.2437;%1.7458;
spwr=sqrt(pwr);
r=1.025;

Ux=1;
Uy=0;
Uz=0;
Uinf=sqrt(Ux*Ux+Uy*Uy+Uz*Uz);
I=5;
up2= (I/100*Uinf).^2;
Rij(1,1) = up2;
Rij(2,2) = up2;
Rij(3,3) = up2;
TKE=3/2*up2;

a11=sqrt(Rij(1,1));
a21=Rij(2,1)/a11;
a22=sqrt( Rij(2,2)-a21^2 );
a31=Rij(3,1)/a11;
a32=(Rij(3,2)-a21*a31)/a22;
a33=sqrt( Rij(3,3)-a31*a31 - a32*a32);


Mx=5;
My=5;
Mz=121;

Dx=10;
Dy=10;
Dz=10;

Lx=1;
Ly=1;
Lz=1;


mx=(Mx-1)/2+1;
my=(My-1)/2+1;
mz=(Mz-1)/2+1;mz=Mz-1;

dx=Dx/(Mx-1); %x grid spacing
dy=Dy/(My-1); %y grid spacing
dz=Dz/(Mz-1); %z grid spacing
dt=dx/Ux; %taylor hypothesis time step

%

dr(1)=1;
for i=2:Mz-1
dr(i)=dr(i-1)*r;
end
R=[0,cumsum(dr)];

y=linspace(0,Dy,My);
drz=dr/dr(mz)*dz; z=[0,cumsum(drz)];
Lskew_factor=sqrt( sum( (drz./dz).^2 )/numel(drz) );

[ymesh,zmesh]=meshgrid(y,z);ymesh=ymesh'; zmesh=zmesh';

mindy=min(y(2:end)-y(1:end-1));
mindz=min(z(2:end)-z(1:end-1));
mzcart=ceil((z(mz)-z(1))/mindz);
mycart=ceil((y(my)-y(1))/mindy);

dycart=(y(my)-y(1))/(mycart-1);
ycart=y(1):dycart:y(end);
dzcart=(z(mz)-z(1))/(mzcart-1);
zcart=[z(1):dzcart:mzcart*dzcart,(mzcart+1)*dzcart:dzcart:z(end)];
[ymeshCart,zmeshCart]=meshgrid(ycart,zcart);ymeshCart=ymeshCart'; zmeshCart=zmeshCart';
Mycart=numel(ycart);
Mzcart=numel(zcart);

ny=ceil(Ly/dy); %number of nodes to sum correlation function over
nz=ceil(Lz/dz);
nx=ceil(Lx/dt);

Nx=2*nx;
Ny=2*ny;
Nz=2*nz;

rng(1);
Rx = sqrt(3)*(2*rand(2*Nx+1,Ny+My+Ny+1, Nz+Mz+Nz+1)-1);
Ry = sqrt(3)*(2*rand(2*Nx+1,Ny+My+Ny+1, Nz+Mz+Nz+1)-1);
Rz = sqrt(3)*(2*rand(2*Nx+1,Ny+My+Ny+1, Nz+Mz+Nz+1)-1);

%filter coefficients
btil = @(k,n) exp( -pi*k.^2/(2*n^2) );
bk = @(k,n,N) btil(k,n)./sqrt( sum( btil(-N:N,n).^2 ) );

btilz = @(k,n) exp( -pi/2*abs(k/n).^pwr );
bkz = @(k,n,N) btilz(k,n)./sqrt( sum( btil(-N:N,n).^2 ) );

btil2 = @(k,n) exp( -pi/2*(r.^(2*k)-2*r.^k+1)/(r^(2*n)-2*r^n+1) );
bk2 = @(k,n,Nm,Np) btil2(k,n)./sqrt( sum( btil(-Nm:Np,n).^2 ) );

%Generate filter coefficients
for i=1:2*Nx+1
  ii=i-Nx-1;
  BI(i) =  bk(ii,nx,Nx);
end

for j=1:2*Ny+1
  jj=j-Ny-1;
  BJ(j) =  bk(jj,ny,Ny);
end


for k=1:2*Nz+1
  kk=k-Nz-1;
  BK(k) =  bkz(kk,nz,Nz);
%   BK(k) =  bk2(kk,nz,Nz,Nz);
end

bi=repmat(BI',[1,Ny+My+Ny+1, Nz+Mz+Nz+1]);
for k=1:2*Nz+1
    for j=1:2*Ny+1
        bjk(j,k) = BJ(j)*BK(k);
        
        for i=1:2*Nx+1
          bijk(i,j,k)=BI(i)*BJ(j)*BK(k);
        end
    end
end


ip=[-Nx:Nx]+(Nx+1);
jp=[-Ny:Ny]+(Ny+1);
kp=[-Nz:Nz]+(Nz+1);


% uu=zeros(Mycart,Mzcart);
% Mycart=My; Mzcart=Mz; mycart=my; mzcart=mz;
% ycart=y; zcart=z;
uu=zeros(Mycart,Mzcart);
vv=uu;
ww=vv;
uv=uu;
vw=uu;

UU_cy = zeros(1, 2*(Mycart-1)+1);
UU_cz = zeros(1, 2*(Mzcart-1)+1);
UU_y = UU_cy;
UU_z = UU_cz;
UU_cza = zeros(1, 2*(Mzcart-1)+1);
UU_za = UU_cza;

cnt=1;

figure(1)


for nn=1:niter


RRx = squeeze(sum(bi.*Rx,1));
RRy = squeeze(sum(bi.*Ry,1));
RRz = squeeze(sum(bi.*Rz,1));
for k=1:Mz
    for j=1:My
        ux(j,k) = sum(sum(bjk.*RRx(j+jp,k+kp)));
        uy(j,k) = sum(sum(bjk.*RRy(j+jp,k+kp)));
        uz(j,k) = sum(sum(bjk.*RRz(j+jp,k+kp)));
    end
end

  Rx(1:end-1,:,:)=Rx(2:end,:,:);
  Ry(1:end-1,:,:)=Ry(2:end,:,:);
  Rz(1:end-1,:,:)=Rz(2:end,:,:);
      
  Rx(end,:,:)= sqrt(3)*(2*rand(2*Ny+My+1,2*Nz+Mz+1)-1);
  Ry(end,:,:)= sqrt(3)*(2*rand(2*Ny+My+1,2*Nz+Mz+1)-1);
  Rz(end,:,:)= sqrt(3)*(2*rand(2*Ny+My+1,2*Nz+Mz+1)-1);

  xVel=Ux+a11*ux+a21*uy+a31*uz;
  yVel=Uy+a22*uy+a32*uz;
  zVel=Uz+a33*uz;

  xVelO=xVel;
  yVelO=yVel;
  zVelO=zVel;
  xVel=interp2(ymesh',zmesh',xVel',ymeshCart',zmeshCart','spline'); xVel=xVel';
  yVel=interp2(ymesh',zmesh',yVel',ymeshCart',zmeshCart','spline'); yVel=yVel';
  zVel=interp2(ymesh',zmesh',zVel',ymeshCart',zmeshCart','spline'); zVel=zVel';
  
  % Evaluate second-order statistics ======================================
  up=xVel-Ux;
  vp=yVel-Uy;
  wp=zVel-Uz;
  
  uu = uu + up.*up;
  vv = vv + vp.*vp;
  ww = ww + wp.*wp;
  uv = uv + up.*vp;
  vw = vw + vp.*wp;
  
  if (mod(nn,50)==0)
      UU(:,:,cnt) = uu/nn;
      VV(:,:,cnt) = vv/nn;
      WW(:,:,cnt) = ww/nn;
      
      UV(:,:,cnt) = uv/nn;
      VW(:,:,cnt) = vw/nn;
      NN(cnt)=nn;
      cnt=cnt+1;
  end

  % Evaluate correlation function =========================================


  for jj=1:Mycart
    j=[jj:Mycart]-mycart;
    uu_y1(jj) = sum(up(:,mzcart).*up(:,mzcart));
    uu_cy1(jj) = sum(up(mycart+j-jj+1,mzcart).*up(mycart+j,mzcart));
    
    uu_cy1(jj) = sum(up(jj:Mycart,mzcart).*up(1:Mycart-jj+1,mzcart));
    
    j=[Mycart:-1:jj]-mycart;
    uu_y2(jj) = sum(up(:,mzcart).*up(:,mzcart));
%     uu_cy2(jj) = sum(up(mycart-j+jj-1,mzcart).*up(mycart-j,mzcart));
    
    uu_cy2(jj) = sum(up(jj:Mycart,mzcart).*up(1:Mycart-jj+1,mzcart));

  end
  for jj=1:Mzcart
    j=[jj:Mzcart]-mzcart;
    uu_z1(jj) = sum(up(mycart,:).*up(mycart,:));
    uu_cz1(jj) = sum(up(mycart,mzcart+j-jj+1).*up(mycart,mzcart+j));
    
%     j=[Mzcart:-1:jj]-mzcart;
    uu_z2(jj) = sum(up(mycart,:).*up(mycart,:));
    uu_cz2(jj) = sum(up(mycart,jj:Mzcart).*up(mycart,1:Mzcart-jj+1));
    
    nnn=numel(jj:Mzcart);
%     plot([1:nnn],up(mycart,mzcart+j-jj+1),'r-o',[1:nnn],up(mycart,mzcart+j),'g-*',[1:nnn],up(mycart,jj:Mzcart).*up(mycart,1:Mzcart-jj+1),'k-o')
    
    zuu(jj) = zmeshCart(mycart,jj)-zmeshCart(mycart,1);
  end
  
  mzcart2 = mzcart;
  mycart2 = mycart-1;
  for jj=1:Mzcart
    j=[jj:Mzcart]-mzcart2;
    uu_z1a(jj) = sum(up(mycart2,:).*up(mycart2,:));
    uu_cz1a(jj) = sum(up(mycart2,mzcart2+j-jj+1).*up(mycart2,mzcart2+j));
    
%     j=[Mzcart:-1:jj]-mzcart2cart;
    uu_z2a(jj) = sum(up(mycart2,:).*up(mycart2,:));
    uu_cz2a(jj) = sum(up(mycart2,jj:Mzcart).*up(mycart2,1:Mzcart-jj+1));
    
    nnn=numel(jj:Mzcart);
%     plot([1:nnn],up(mycart2,mzcart2+j-jj+1),'r-o',[1:nnn],up(mycart2,mzcart2+j),'g-*',[1:nnn],up(mycart2,jj:Mzcart).*up(mycart2,1:Mzcart-jj+1),'k-o')
    
    zuua(jj) = zmeshCart(mycart2,jj)-zmeshCart(mycart2,1);
  end
  

%   for jj=1:Mx
%       j=[jj:Mx]-mx;
%     uu_1p(jj) = sum(up(:,my).*up(:,my));
%     uu_cxp(jj) = sum(up(mx+j-jj+1,my).*up(mx+j,my));
%     uu_2p(jj) = sum(up(mx,:).*up(mx,:));
%     uu_cyp(jj) = sum(up(mx,my+j-jj+1).*up(mx,my+j));
%     
%     j=[Mx:-1:jj]-mx;
%     uu_1m(jj) = sum(up(:,my).*up(:,my));
%     uu_cxm(jj) = sum(up(mx-j+jj-1,my).*up(mx-j,my));
%     uu_2m(jj) = sum(up(mx,:).*up(mx,:));
%     uu_cym(jj) = sum(up(mx,my-j+jj-1).*up(mx,my-j));
%   end
  uu_cy = [fliplr(uu_cy2(2:end)),uu_cy1];
  uu_cz = [fliplr(uu_cz2(2:end)),uu_cz1];
  uu_y = [fliplr(uu_y2(2:end)),uu_y1];
  uu_z = [fliplr(uu_z2(2:end)),uu_z1];
  
  uu_cza = [fliplr(uu_cz2a(2:end)),uu_cz1a];
  uu_za = [fliplr(uu_z2a(2:end)),uu_z1a];
  
  UU_cy = UU_cy + uu_cy;
  UU_cz = UU_cz + uu_cz;
  UU_y = UU_y + uu_y;
  UU_z = UU_z + uu_z;
  
  UU_cza = UU_cza + uu_cza;
  UU_za = UU_za + uu_za;
  
  Ruu=@(x,L,x0) exp(-pi*(x-x0).^pwr/(4*(Lskew_factor*L).^pwr));
  plot(zuu,Ruu(zmeshCart(mycart,:),Lz,0),'k--',...
       zuu, uu_cz(end-Mzcart+1:end)./uu_z(end-Mzcart+1:end),'r-o',...
       zuu, UU_cz(end-Mzcart+1:end)./UU_z(end-Mzcart+1:end),'g-*',...
       zuua,uu_cz1a./uu_z1a,'b-v',...
       zuua,UU_cza(end-Mzcart+1:end)./UU_za(end-Mzcart+1:end),'c-^');
  pause(0.001)
  

  
%   if (mod(nn,100)==0)
% %       figure(1)
% %       subplot(2,3,1)
% %       surf(ymeshCart',zmeshCart',xVel)
% % %       surf(uu/nn)
% %       view([0,0,90])
% %       shading interp
% %       axis equal
% %       title('x velocity')
% %       
% %       subplot(2,3,4)
% %       surf(ymesh',zmesh',xVelO)
% % %       surf(uu/nn)
% %       view([0,0,90])
% %       shading interp
% %       axis equal
% %       title('x velocity')
% % 
% % 
% %       subplot(2,3,2)
% %       surf(ymeshCart',zmeshCart',yVel)
% % %       surf(vv/nn)
% %       view([0,0,90])
% %       shading interp
% %       axis equal
% %       title('y velocity')
% % 
% %       subplot(2,3,5)
% %       surf(ymesh',zmesh',yVelO)
% % %       surf(vv/nn)
% %       view([0,0,90])
% %       shading interp
% %       axis equal
% %       title('y velocity')
% %       
% %       subplot(2,3,3)
% %       surf(ymeshCart',zmeshCart',zVel)
% % %       surf(ww/nn)
% %       view([0,0,90])
% %       shading interp
% %       axis equal
% %       title('z velocity')
% %       
% %       subplot(2,3,6)
% %       surf(ymesh',zmesh',zVelO)
% % %       surf(ww/nn)
% %       view([0,0,90])
% %       shading interp
% %       axis equal
% %       title('z velocity')
% 
% 
%     UU_cy1=UU_cy/nn;
%     UU_cz1=UU_cz/nn;
%     UU_y1=UU_y/nn;
%     UU_z1=UU_z/nn;
% 
%     jj=(1:numel(y));
%     kk=(1:numel(z));
%     
%     yycart=[fliplr(-ycart(2:end)),ycart];
%     zzcart=[fliplr(-zcart(2:end)),zcart];
%     
%     ny_new=floor(Ly/dycart);
%     nz_new=floor(Lz/dzcart);
%     Ruud=@(k,n,i) exp(-pi*(k-i).^2/(4*n.^2));
%     Ruu=@(x,L,x0) exp(-pi*(x-x0).^2/(4*L.^2));
%     errDy=max(abs((Ruud(jj,ny,my)-Ruu(y,Ly,y(my)))))*100;
%     errDz=max(abs((Ruud(kk,nz,mz)-Ruu(z,Lz,z(mz)))))*100;
%     err_Ruuy = max(abs((UU_cy1./UU_y1-Ruu(yycart,Ly,0))))*100;
%     err_Ruuz = max(abs((UU_cz1./UU_z1-Ruu(zzcart,Lz,0))))*100;
%     
%     figure(2)
%     subplot(2,1,1)
% 
%     
%     plot(yycart+ycart(mycart),UU_cy1./UU_y1, 'k-o',...
%          ycart, Ruu(ycart,Ly,ycart(mycart)),'k--',y, Ruud(jj,ny,my),'b-.*')
%     text(ycart(mycart)+3,0.8,['Error: Ruu disc: ',num2str(errDy,5)])
%     text(ycart(mycart)+3,0.7,['Error: Ruu y: ',num2str(err_Ruuy,5)])
%     
%     subplot(2,1,2)
%     plot(zzcart+zcart(mzcart),UU_cz1./UU_z1, 'k-o',...
%         zcart, Ruu(zcart,Lz,zcart(mzcart)),'k--',z, Ruud(kk,nz,mz),'b-.*')
%     text(zcart(mzcart)+3,0.8,['Error: Ruu disc: ',num2str(errDz,5)])
%     text(zcart(mzcart)+3,0.7,['Error: Ruu z: ',num2str(err_Ruuz,5)])
%     
%     pause(0.1)
%   
% 
%   end
  
end

    Ruu_out = UU_cz(end-Mzcart+1:end)./UU_z(end-Mzcart+1:end);
     Ruu_errfun = @(n) max( abs( (exp(-pi/4*(zuu/(Lskew_factor*Lz)).^n)-Ruu_out) ) ); 
     pwr_out=fminsearch(Ruu_errfun,pwr);
%      nvec=linspace(1,2,20);
%      for j=1:numel(nvec)
%          errfun(j) = Ruu_errfun(nvec(j));         
%      end
%      figure
      plot(zuu,exp(-pi/4*(zuu/(Lskew_factor*Lz)).^pwr_out),'k-o', zuu, Ruu_out,'r-o')

puinstauu=uu/niter;
vv=vv/niter;
ww=ww/niter;

k=1/2*(uu+vv+ww);
Intensity=sqrt(2/3*k);

for j=1:numel(NN)
    UUm(j)=mean(mean(UU(:,:,j)));
    VVm(j)=mean(mean(VV(:,:,j)));
    WWm(j)=mean(mean(WW(:,:,j)));
    UVm(j)=mean(mean(UV(:,:,j)));
    VWm(j)=mean(mean(VW(:,:,j)));
end

UU_cy1=UU_cy/niter;
UU_cz1=UU_cz/niter;
UU_y1=UU_y/niter;
UU_z1=UU_z/niter;

yycart=[fliplr(-ycart(2:end)),ycart];
zzcart=[fliplr(-zcart(2:end)),zcart];

ny_new=floor(Ly/dycart);
nz_new=floor(Lz/dzcart);
Ruud=@(k,n,i) exp(-pi*(k-i).^2/(4*n.^2));
Ruu=@(x,L,x0) exp(-pi*(x-x0).^2/(4*L.^2));
errDy=max(abs((Ruud(jj,ny,my)-Ruu(y,Ly,y(my)))))*100;
errDz=max(abs((Ruud(kk,nz,mz)-Ruu(z,Lz,z(mz)))))*100;
err_Ruuy = max(abs((UU_cy1./UU_y1-Ruu(yycart,Ly,0))))*100;
err_Ruuz = max(abs((UU_cz1./UU_z1-Ruu(zzcart,Lz,0))))*100;

figure(2)
subplot(2,1,1)


plot(yycart+ycart(mycart),UU_cy1./UU_y1, 'k-o',...
     ycart, Ruu(ycart,Ly,ycart(mycart)),'k--',y, Ruud(jj,ny,my),'b-.*')
text(ycart(mycart)+3,0.8,['Error: Ruu disc: ',num2str(errDy,5)])
text(ycart(mycart)+3,0.7,['Error: Ruu y: ',num2str(err_Ruuy,5)])

subplot(2,1,2)
plot(zzcart+zcart(mzcart),UU_cz1./UU_z1, 'k-o',...
    zcart, Ruu(zcart,Lz,zcart(mzcart)),'k--',z, Ruud(kk,nz,mz),'b-.*')
text(zcart(mzcart)+3,0.8,['Error: Ruu disc: ',num2str(errDz,5)])
text(zcart(mzcart)+3,0.7,['Error: Ruu z: ',num2str(err_Ruuz,5)])