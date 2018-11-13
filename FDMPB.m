function sol_u=FDMPB(V,N,eps)
V_L=-V;
V_R=V;
deltax=2/N;
lambda=eps^2/(deltax^2);
x=-1:deltax:1;
phi_old=x*V_R;% phi_old also satisfys the boundary conditions
w=zeros(N-1,1);
phi=zeros(N+1,1);
error3=1;
phi(1)= V_L;% boundary conditions
phi(N+1) =V_R;
B0 =sparse([1:1:N-1,1:1:N-2 ,2:1:N-1],[1:1:N-1, 2:1:N-1, 1:1:N-2],...
    [2*ones(1,N-1),-ones(1,N-2) ,-ones(1,N-2)],N-1,N-1)*lambda;
while error3>10^(-11)
    B =B0+sparse(1:N-1,1:N-1,cosh(phi_old(2:N)),N-1,N-1);
    w(1:N-1)=cosh(phi_old(2:N)).*phi_old(2:N)-sinh(phi_old(2:N));
    w(1)=w(1)+lambda*phi(1);% use the boundary condition
    w(N-1)=w(N-1)+lambda*phi(N+1);
    phi(2:N)=B\w;
    EError3=phi-phi_old';
    error3=max(abs(EError3));
    phi_old=phi';
end
sol_u=phi;