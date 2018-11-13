function [u0,timerbm]=cap_rbm(Nx,epsc)
load al400
tic
hx = 2/Nx;
Ny=Nx;
hy = 2/Ny;
xx=(-1+hx:hx:1-hx);
TN=(Ny+1)*(Nx-1); % unknown nodes
for i3 = 1:Ny+1
    F1(1,(i3-1)*(Nx-1)+1:(i3-1)*(Nx-1)+Nx-1)=gaussmf(xx,[0.1 0])*gaussmf(-1+hy*(i3-1),[0.1 0]);
end
qq=zeros(102,Ny);
Cap_rbm=zeros(101,length(epsc));
PHI=zeros(Nx+1,Ny+1);
for ic=1:length(epsc)
    eps_test=epsc(ic);
    lambdax=(eps_test/hx)^2;
    lambday =(eps_test/hy)^2;
    A1 = ALL1*lambdax;
    A2 = ALL2*lambday;
    A12 = A1+A2;
    
    for i1= 0:101
        v_test= (i1-1)*0.02;
        for i3 = 1:Ny+1
            uold_test(1,(i3-1)*(Nx-1)+1:(i3-1)*(Nx-1)+Nx-1) = v_test*xx;
        end
        V_L = - v_test;
        V_R = v_test;
        w=zeros(TN,1);
        errPB=1;
        itertest = 1;
        while errPB>10^(-8)
            %% 这个AD表示主对角线上那一部分产生的结果，其实只有主对角线上的元素
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
            itertest=itertest+1;
        end
        PHI(1,1:Ny+1)= V_L;% boundary condition
        PHI(Nx+1,1:Ny+1)= V_R;% boundary condition
        for j0=1:Ny+1
            PHI(2:Nx,j0)=u_test((j0-1)*(Nx-1)+1:Nx-1+(j0-1)*(Nx-1));
        end
        qq(i1+1,:)=eps_test^2*(4*PHI(2,1:Ny)-3*PHI(1,1:Ny)-PHI(3,1:Ny))*hy*Nx/4;
%         qq_fdm(i1+1,:)=capFDM(Nx,eps_test,v_test,F1);
    end
    Q=sum(qq,2);
%     Q_fdm=sum(qq_fdm,2);
    Cap_rbm(:,ic)=(Q(2:102)-Q(1:101))/0.02/2;
%     Cap_fdm(:,ic)=(Q_fdm(2:102)-Q_fdm(1:101))/0.02/2;
end
timerbm=toc;
u0=Cap_rbm;
% Err=Cap_rbm-Cap_fdm;
% err_abs=abs(Err);

%%


