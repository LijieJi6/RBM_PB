% capacitance
clear
load al200
tic
Nx=200;
Ny=Nx;
preN=20;
hx = 2/Nx;
hy = 2/Ny;
xx=(-1+hx:hx:1-hx);%不含有边界，边界电压已知
TN=(Ny+1)*(Nx-1);%unknown nodes
TN1=(Ny-1)*(Nx-1);
eps_test=0.25;
v_test=3.57;

V_L = - v_test;
V_R = v_test;
lambdax=(eps_test/hx)^2;
lambday =(eps_test/hy)^2;
w=zeros(TN,1);
A1 = ALL1*lambdax;
A2 = ALL2*lambday;
A12 = A1+A2;
for i3 = 1:Ny+1
    uold_test(1,(i3-1)*(Nx-1)+1:(i3-1)*(Nx-1)+Nx-1) = v_test*xx;
    F1(1,(i3-1)*(Nx-1)+1:(i3-1)*(Nx-1)+Nx-1)=gaussmf(xx,[0.1 0])*gaussmf(-1+hy*(i3-1),[0.1 0]);
end
errPB=1;
while errPB>10^(-8)
    BB = sparse((1:1:TN),(1:1:TN),...
        cosh(uold_test),TN,TN);
    AD=temxi*(BB*temxi');
    B = A12 +AD;
    %% define right hand vector
    w(:,1)=cosh(uold_test).*uold_test-sinh(uold_test)+F1;
    w(1:(Nx-1):(Nx-1)*Ny+1,1)=w(1:(Nx-1):(Nx-1)*Ny+1,1)+lambdax*V_L;
    w(Nx-1:(Nx-1):TN)=w(Nx-1:(Nx-1):TN,1)+lambdax*V_R;
    F=temxi*w;
    coeff_test  = B\F;
    % store all the nodes information of all paras
    u_test=coeff_test'*temxi;
    errPB=max(abs((u_test-uold_test)'));
    uold_test=u_test;
end
% BC=[(1:Nx+1:1+Ny*(Nx+1)),(Nx+1:Nx+1:(Ny+1)*(Nx+1))];
% RBM_u(setdiff((1:1:(Ny+1)*(Nx+1)),BC))=u_test;
% RBM_u(1:Nx+1:1+Ny*(Nx+1))= V_L;% boundary condition
% RBM_u(Nx+1:Nx+1:(Ny+1)*(Nx+1))= V_R;% boundary condition

toc

%{

%RBM_u是所有节点矩阵值，电势在X方向为V_L到V_R的PB电势曲线
%也即每一个Ny节点对应的列向量是一条PB电势曲线
%以下代码是为了方便画出来我们的电势图，mesh


PHI(1,1:Ny+1)= V_L;% boundary condition
PHI(Nx+1,1:Ny+1)= V_R;% boundary condition
for j0=1:Ny+1
    PHI(2:Nx,j0)=u_test((j0-1)*(Nx-1)+1:Nx-1+(j0-1)*(Nx-1));
end
% mesh(-1:hy:+1,-1:hx:+1,(PHI),'LineWidth',3)
mesh(-1:hy:+1,-1:hx:+1,PHI,'LineWidth',3)
pos=axis;%[xmin xmax ymin ymax zmin zmax]
% xlabel('\fontsize{30} \it y','FontName', 'Times New Roman','position',[0.1,-1.1,0])
% ylabel('\fontsize{30} \it x','FontName', 'Times New Roman','position',[-1.1,0.2,0])
% zlabel('\fontsize{19} \it Potential Distribution','FontName', 'Times New Roman','position',[-1.35,1,0.15])

% xlabel('\fontsize{30} \it y','FontName', 'Times New Roman','position',[0.1,-1.1,-0.2])
% ylabel('\fontsize{30} \it x','FontName', 'Times New Roman','position',[-1.1,0.25,-0.2])
% zlabel('\fontsize{19} \it Potential Distribution','FontName', 'Times New Roman','position',[-1.4,1,0.15])

% xlabel('\fontsize{30} \it y','FontName', 'Times New Roman','position',[0.1,-1.1,-0.5])
% ylabel('\fontsize{30} \it x','FontName', 'Times New Roman','position',[-1.1,0.25,-0.5])
% zlabel('\fontsize{18} \it Potential Distribution','FontName', 'Times New Roman','position',[-1.4,1,0.05])

xlabel('\fontsize{30} \it y','FontName', 'Times New Roman','position',[0.1,-1.1,-5])
ylabel('\fontsize{30} \it x','FontName', 'Times New Roman','position',[-1.1,0.25,-5])
zlabel('\fontsize{19} \it Potential Distribution','FontName', 'Times New Roman','position',[-1.33,1,0.1])
%title('Potential distribution in [-1,1] \times [-1,1]')
set(gca,'FontName','Times New Roman','FontSize',18,...
    'GridColor','k','FontWeight','bold','LineWidth',2)
%set(gca,'position',[0.2 0.2 0.5 0.5]);
%axis([-1 1 -1 1 -0.5 0.5])
box on
colormap(cool)
% colormap(rand(64,3))
% 0  0.1 1 5


toc
tic
u_accurate=FDMtwo(Nx,Ny,eps_test,v_test);% whole solution
toc
max(abs(u_accurate-RBM_u))

%}
