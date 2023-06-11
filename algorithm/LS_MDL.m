function [Signaldim,MDL] = LS_MDL(D,N,L)
D = sort(D,'descend');
m = 30;
% m = min(N,L) - 1;
MDL = zeros(1, m - 1);
for k = 0:m - 1
    tau = 1/(m - k)*sum(D(k+1:m));
    beta = (sum(D(k+1:m).^2) + sum(D(k+1:m))^2)/((N+1)*(sum(D(k+1:m).^2) - sum(D(k+1:m))^2/(m-k)));
    beta = min(beta, 1);
    rho = zeros(1,m - k);
    for i = k+1:m
        rho(i - k) = beta*tau + (1 - beta)*D(i);
    end
    MDL(k + 1) = N*(m-k)*log(1/(m-k)*sum(rho)/nthroot(prod(rho),m-k)) + 1/2*k*(k-1)*log(N);
end
[~,Signaldim] = min(MDL);