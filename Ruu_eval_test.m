
clc
clear all
L=2;


Ruu=@(x,L) exp(-pi*x.^2/(4*L.^2));
b=@(x,L) exp(-pi*x.^2/(2*L.^2));

dx=0.1;
N=4;
nn=2*N/dx;
x=-N*L:dx:N*L;

f1=b(x,L);

uu=sum(f1.*f1);


for i=1:numel(x)
    x0 = dx*(i-nn-1);
    f2(i) = sum(  f1.*b(x+x0,L) );
    f3(i) = sum(  f1(i:end).*f1(1:end-i+1) );
    
end


out = f2./uu;
out3 = f3./uu;
plot(x,out,'k-o',x+2*N,out3,'g-*',x,Ruu(x,L),'r--');


