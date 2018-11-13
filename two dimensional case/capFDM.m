function u0=capFDM(Nx,eps_test,v_test,F1)
Ny=Nx;
hx=2/Nx;
hy=2/Ny;%这里先设置等分相等
xx=-1+hx:hx:1-hx;
lambdax=(eps_test/hx)^2;
lambday=(eps_test/hy)^2;
TN=(Ny+1)*(Nx-1);%total rows TN
PHI=zeros(Nx+1,Ny+1);
F=zeros(TN,1);
%%
secd1=(1:1:TN-1);
secd2=(Nx-1:Nx-1:Ny*(Nx-1));
secd3= Nx:1:Ny*(Nx-1);
secd4 = 1:1:Nx-1;% modify 3 4
SECD12 = -lambdax*ones(1,TN-1-Ny);% matrix element
SECD34 = -lambday*ones(1,(Ny-1)*(Nx-1));
SECDM34 = -2*lambday*ones(1,Nx-1);
%%
Coeff0 = sparse([1:1:TN ,setdiff(secd1,secd2),setdiff(secd1+1,secd2+1)...
    ,secd3,secd3,secd4,(Nx-1)*Ny+secd4],...
    [1:1:TN ,setdiff(secd1+1,secd2+1),setdiff(secd1,secd2),...
    secd3+Nx-1,secd3-Nx+1,secd4+Nx-1,(Nx-1)*(Ny-1)+secd4],...
    [(2*lambdax+2*lambday)*ones(1,TN),SECD12,SECD12,...
    SECD34,SECD34,SECDM34, SECDM34],...
    TN,TN);

V= v_test;
V_L=-V;
V_R=V;
for iy=1:Ny+1
    phi_old(1,(iy-1)*(Nx-1)+1:(iy-1)*(Nx-1)+Nx-1)=V_R*xx;
end
error=1;
iter=0;
while error>10^(-11)
    Coeff1 = Coeff0+ sparse(1:1:TN ,1:1:TN ,cosh(phi_old),TN,TN);
    F(:,1)=cosh(phi_old).*phi_old-sinh(phi_old)+F1;
    F(1:(Nx-1):(Nx-1)*Ny+1,1)=F(1:(Nx-1):(Nx-1)*Ny+1,1)+lambdax*V_L;
    F(Nx-1:(Nx-1):TN)=F(Nx-1:(Nx-1):TN,1)+lambdax*V_R;
    %Coeff1*phi=F
    phi=Coeff1\F;
    error=max(abs(phi'-phi_old));
    phi_old=phi';
    iter=iter+1;
end
PHI(1,1:Ny+1)= V_L;% boundary condition
PHI(Nx+1,1:Ny+1)= V_R;% boundary condition
for j0=1:Ny+1
    PHI(2:Nx,j0)=phi((j0-1)*(Nx-1)+1:Nx-1+(j0-1)*(Nx-1));
end
u0=eps_test^2*(4*PHI(2,1:Ny)-3*PHI(1,1:Ny)-PHI(3,1:Ny))*hy*Nx/4;


