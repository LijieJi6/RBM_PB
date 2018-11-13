
clear
tic
V_test=3.75;
eps_test=0.2;
N=1000;
%%
deltax=2/N;
lambda=eps_test^2/(deltax^2);
x=-1:deltax:1;
%%
phi_old=x*V_test;% phi_old also satisfys the boundary conditions
w=zeros(N-1,1);
phi=zeros(N+1,1);
error3=1;
phi(1)= -V_test;% boundary conditions
phi(N+1) =V_test;
B0 =sparse([1:1:N-1,1:1:N-2 ,2:1:N-1],[1:1:N-1, 2:1:N-1, 1:1:N-2],...
    [2*ones(1,N-1),-ones(1,N-2) ,-ones(1,N-2)],N-1,N-1)*lambda;
while error3>10^(-11)
    B =B0+sparse(1:N-1,1:N-1,cosh(phi_old(2:N)),N-1,N-1);
    w(1:N-1)=cosh(phi_old(2:N)).*phi_old(2:N)-sinh(phi_old(2:N));
    w(1)=w(1)+lambda*phi(1);% use the boundary condition
    w(N-1)=w(N-1)+lambda*phi(N+1);
    phi(2:N)=B\w;
    error3=max(abs(phi-phi_old'));
    phi_old=phi';
end


toc