%DNS BC Test results
% Plots data from executable Test created by DNSbcTest.f90
clear all 
clc


M1=importdata('VelocityTest.txt');
M2=importdata('SecondMomentsTest.txt');
uL=importdata('Xvelocity_lengthscale.txt');
vL=importdata('Yvelocity_lengthscale.txt');
wL=importdata('Zvelocity_lengthscale.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Average Velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(3,1,1)
plot(M1(:,1),M1(:,2),'r-o',M1(:,1),M1(:,3),'k-*')
title('X velocity')

subplot(3,1,2)
plot(M1(:,1),M1(:,4),'r-o',M1(:,1),M1(:,5),'k-*')
title('Y velocity')

subplot(3,1,3)
plot(M1(:,1),M1(:,6),'r-o',M1(:,1),M1(:,7),'k-*')
title('Z velocity')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Perturbation Velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
subplot(3,2,1)
plot(M2(:,1),M2(:,2),'r-o',M2(:,1),M2(:,3),'k-*')
title('uu')

subplot(3,2,3)
plot(M2(:,1),M2(:,4),'r-o',M2(:,1),M2(:,5),'k-*')
title('vv')

subplot(3,2,5)
plot(M2(:,1),M2(:,6),'r-o',M2(:,1),M2(:,7),'k-*')
title('ww')

subplot(3,2,2)
plot(M2(:,1),M2(:,8),'r-o',M2(:,1),M2(:,9),'k-*')
title('uv')

subplot(3,2,4)
plot(M2(:,1),M2(:,10),'r-o',M2(:,1),M2(:,11),'k-*')
title('uw')

subplot(3,2,6)
plot(M2(:,1),M2(:,12),'r-o',M2(:,1),M2(:,13),'k-*')
title('vw')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Length Scales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
subplot(3,3,1)
plot(uL(:,1),uL(:,2),'r-o',uL(:,1),uL(:,3),'k-*')
title('X-velocity: Lx')

subplot(3,3,2)
nn=numel(uL(:,5));
l2=sqrt( sum( uL(:,5).^2 )/nn );
plot(uL(:,1),uL(:,4),'r-o',uL(:,1),uL(:,5),'k-*', uL(:,1), l2*ones*size(nn,1),'k-v')
title('X-Velocity: Ly')

subplot(3,3,3)
plot(uL(:,1),uL(:,6),'r-o',uL(:,1),uL(:,7),'k-*')
title('X-Velocity: Lz')

%%%%%%%%%%%%%%
subplot(3,3,4)
plot(vL(:,1),vL(:,2),'r-o',vL(:,1),vL(:,3),'k-*')
title('Y-velocity: Lx')

subplot(3,3,5)
nn=numel(vL(:,5));
l2=sqrt( sum( vL(:,5).^2 )/nn );
plot(vL(:,1),vL(:,4),'r-o',vL(:,1),vL(:,5),'k-*', vL(:,1), l2*ones*size(nn,1),'k-v')
title('Y-Velocity: Ly')

subplot(3,3,6)
plot(vL(:,1),vL(:,6),'r-o',vL(:,1),vL(:,7),'k-*')
title('Y-Velocity: Lz')

%%%%%%%%%%%%%%
subplot(3,3,7)
plot(wL(:,1),wL(:,2),'r-o',wL(:,1),wL(:,3),'k-*')
title('Z-velocity: Lx')

subplot(3,3,8)
nn=numel(wL(:,5));
l2=sqrt( sum( wL(:,5).^2 )/nn );
plot(wL(:,1),wL(:,4),'r-o',wL(:,1),wL(:,5),'k-*', wL(:,1), l2*ones*size(nn,1),'k-v')
title('Z-Velocity: Ly')

subplot(3,3,9)
plot(wL(:,1),wL(:,6),'r-o',wL(:,1),wL(:,7),'k-*')
title('Z-Velocity: Lz')