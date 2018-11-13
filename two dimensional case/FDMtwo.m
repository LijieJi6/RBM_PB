function sol_u = FDMtwo(Nx,Ny,eps,V)
V_L=-V;
V_R=V;
hx=2/Nx;
hy=2/Ny;
xx=(-1+hx:hx:1-hx);
lambdax=(eps/hx)^2;
lambday=(eps/hy)^2;
for iy=1:Ny+1
    phi_old(1,(iy-1)*(Nx-1)+1:(iy-1)*(Nx-1)+Nx-1)=V_R*xx;
    F1(1,(iy-1)*(Nx-1)+1:(iy-1)*(Nx-1)+Nx-1)=gaussmf(xx,[0.1 0])*gaussmf(-1+hy*(iy-1),[0.1 0]);
end
TN=(Ny+1)*(Nx-1);%total rows TN
TN1=(Ny-1)*(Nx-1);
F=zeros(TN,1);
error=1;
iter=0;
%%
secd1=(1:1:TN-1);
secd2=(Nx-1:Nx-1:Ny*(Nx-1));
secd3= Nx:1:Ny*(Nx-1);
secd4 = 1:1:Nx-1;% modify 3 4
SECD12 = -lambdax*ones(1,TN-1-Ny);% matrix element
SECD34 = -lambday*ones(1,TN1);
SECDM34 = -2*lambday*ones(1,Nx-1);
Coeff0 = sparse([1:1:TN ,setdiff(secd1,secd2),setdiff(secd1+1,secd2+1)...
    ,secd3,secd3,secd4,(Nx-1)*Ny+secd4],...
    [1:1:TN ,setdiff(secd1+1,secd2+1),setdiff(secd1,secd2),...
    secd3+Nx-1,secd3-Nx+1,secd4+Nx-1,(Nx-1)*(Ny-1)+secd4],...
    [(2*lambdax+2*lambday)*ones(1,TN),SECD12,SECD12,...
    SECD34,SECD34,SECDM34, SECDM34],...
    TN,TN);
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
%sol_u=phi'; 
% sol_u in following form is the whole solution
BC=[(1:Nx+1:1+Ny*(Nx+1)),(Nx+1:Nx+1:(Ny+1)*(Nx+1))];
sol_u(setdiff((1:1:(Ny+1)*(Nx+1)),BC))=phi';
sol_u(1:Nx+1:1+Ny*(Nx+1))= V_L;% boundary condition
sol_u(Nx+1:Nx+1:(Ny+1)*(Nx+1))= V_R;% boundary condition


