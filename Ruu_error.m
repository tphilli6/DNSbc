function Ruu_error
%Discrete error approximation for DNS filter coefficients
n=10;

C=1:.1:3;
N=1:20;

for n=1:numel(N)
%     for i=1:numel(C)
%         ME(i) = getmaxerror(N, C(i), 0);
%     end
    
    ME(n) = getmaxerror(N(n), 2, 1);
%     semilogy(N,ME,'k-o',C,0.001*ones(size(C)),'k--')
%     hold on
end
semilogy(N,ME,'k-o',N,0.001*ones(size(N)),'k--')
end


function maxerror = getmaxerror(n,c, ploton)

    Ruufun = @(k,n) exp(-pi/4*k.^2/n.^2);
    bkfun = @(k,n) exp( -pi/2*k.^2/n.^2 );


    N=ceil(c*n);

    k=-N:N;

    bk=bkfun(k,n);
    for i=1:numel(bk)
       BK(i) = sum( bk(1:end+1-i).*bk(i:end) );
    end
    BK=BK/sum( bk.^2 );

    K=(1:numel(BK))-1;
    error = Ruufun(K,n)-BK;
    maxerror = max( abs(error) );

    if ploton
        subplot(2,1,1)
        plot(K,Ruufun(K,n),'k',K,BK,'r-o')
        title('Function Comparison')

        subplot(2,1,2)
        plot(K, Ruufun(K,n)-BK,'r-o')
        title('Error in Function')
        text(1,maxerror*.9,['Max Error = ', num2str(maxerror)])
    end

end