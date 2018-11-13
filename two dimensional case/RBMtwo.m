function u0=RBMtwo(Paraeps,Parav,Nx,Nmax,testeps,testv)
% % offline
Ny=Nx;
hx = 2/Nx;
hy = 2/Ny;
xx=(-1+hx:hx:1-hx);
TN=(Ny+1)*(Nx-1);%unknown nodes
TN1=(Ny-1)*(Nx-1);
lengthx=length(Paraeps);%length of eps
lengthy=length(Parav);%length of potential v
length0=lengthx*lengthy;% total parameter number
preeps=zeros(1,Nmax);
prev=zeros(1,Nmax);

preeps(1)=0.24;
prev(1)=2.5;
for iy=1:Ny+1
    F1(1,(iy-1)*(Nx-1)+1:(iy-1)*(Nx-1)+Nx-1)=gaussmf(xx,[0.1 0])*gaussmf(-1+hy*(iy-1),[0.1 0]);
end
%accurate solution of the first selected parameter
xi(1,:)=FDMtwo(Nx,Ny,preeps(1),prev(1));% whole solution
xi(1,:)=xi(1,:)/sqrt(dot(xi(1,:),xi(1,:)));% normalized basis
BC=[(1:Nx+1:1+Ny*(Nx+1)),(Nx+1:Nx+1:(Ny+1)*(Nx+1))];
temxi(1,:) = xi(1,(setdiff((1:1:(Ny+1)*(Nx+1)),BC)));
%%
eigval=zeros(1,length0);
u_old =zeros(length0,TN);
errRB=zeros(1,Nmax);% eror estimator of the selection
temu=zeros(length0,TN);
F=zeros(TN,1);
temtemu=zeros(1,2);% boundary conditions
%%
secd1=(1:1:TN-1);
secd2=(Nx-1:Nx-1:Ny*(Nx-1));
secd3= Nx:1:Ny*(Nx-1);
secd4 = 1:1:Nx-1;% modify 3 4
%%
CoeffMatrixL1 = sparse([1:1:TN ,setdiff(secd1,secd2),setdiff(secd1+1,secd2+1)],...
    [1:1:TN ,setdiff(secd1+1,secd2+1),setdiff(secd1,secd2)],...
    [2*ones(1,TN),- ones(1,TN-1-Ny),- ones(1,TN-1-Ny)],...
    TN,TN);
ALL1=temxi(1,:)*(CoeffMatrixL1*temxi(1,:)');
CoeffMatrixL2 = sparse([1:1:TN,secd3,secd3,secd4,(Nx-1)*Ny+secd4],...
    [1:1:TN ,2*Nx-1:1:TN,1:1:TN1,secd4+Nx-1,TN1+secd4],...
    [2*ones(1,TN),-ones(1,TN1),-ones(1,TN1),-2*ones(1,Nx-1),-2*ones(1,Nx-1)],...
    TN,TN);
ALL2=temxi(1,:)*(CoeffMatrixL2*temxi(1,:)');
%% initial potential for all parameters
eig0 =eigs((CoeffMatrixL1+CoeffMatrixL2)'*(CoeffMatrixL1+CoeffMatrixL2),1,'sm');

for i1=1:lengthx % select a eps
    for i3 = 1:Ny+1
        u_old((i1-1)*lengthy+1:(i1-1)*lengthy+lengthy,(i3-1)*(Nx-1)+1:(i3-1)*(Nx-1)+Nx-1)...
            =Parav'*xx;
    end
    eigval((i1-1)*lengthy+1:(i1-1)*lengthy+lengthy)=Paraeps(i1)^4/hy^4*eig0;
end
preN=1;
while preN<20
    mf=zeros(length0,TN);
    k0=0;
    for ix=1:lengthx
        lambdax=(Paraeps(ix)/hx)^2;
        lambday =(Paraeps(ix)/hy)^2;
        AL1 =ALL1*lambdax;
        AL2 =ALL2*lambday;
        A12 = AL1+AL2;
        for iy=1:lengthy
            iter=0;
            errPB=1;
            k0=k0+1;
            temtemu(1)= -Parav(iy);
            temtemu(2)= +Parav(iy);% boundary condition
            while errPB>10^(-8)
                BB = sparse((1:1:TN),(1:1:TN),...
                    cosh(u_old(k0,:)),TN,TN);
                AD=temxi(1:preN,:)*(BB*temxi(1:preN,:)');
                B = A12 +AD;
                % right hand vector
                F(:,1)=cosh(u_old(k0,:)).*u_old(k0,:)-sinh(u_old(k0,:))+F1;
                F(1:(Nx-1):(Nx-1)*Ny+1,1)=F(1:(Nx-1):(Nx-1)*Ny+1,1)+lambdax*temtemu(1);
                F(Nx-1:(Nx-1):TN)=F(Nx-1:(Nx-1):TN,1)+lambdax*temtemu(2);
                w=temxi(1:preN,:)*F;
                % solve this linear equation
                coeff=B\w;
                % store all the nodes information of all paras
                temu(k0,:)=coeff'*temxi(1:preN,:);
                errPB=max(abs(temu(k0,:)-u_old(k0,:)));
                u_old(k0,:)=temu(k0,:);
                iter=iter+1;
            end
            % eigenvalues of EB
            EB =lambdax*CoeffMatrixL1+lambday*CoeffMatrixL2+BB;
            mf(k0,:) = EB*temu(k0,:)'-F;% error of each parameter
        end
    end
    delta=sqrt(sum((mf.^2),2))./sqrt(eigval)';% L2 norm£¬£¬£¬
    [x,y]=max(delta);% y is our wanted index, x is the maximum errRB of delta
    %select  a new parameter
    if mod(y,lengthy)~=0
        i2=mod(y,lengthy); i1=floor(y/lengthy)+1;
    else
        i1=floor(y/lengthy); i2=lengthy;
    end
    preeps(preN+1)=Paraeps(i1);
    prev(preN+1)=Parav(i2);
    errRB(preN)=x;
    preN=preN+1
    xi(preN,:)= FDMtwo(Nx, Ny,preeps(preN),prev(preN));
    % GS
    sum0=0;
    for m2=1:preN-1
        sum0=sum0+dot(xi(preN,:),xi(m2,:))/dot(xi(m2,:),xi(m2,:))*xi(m2,:);
    end
    xi(preN,:)=xi(preN,:)-sum0;
    xi(preN,:)=xi(preN,:)/sqrt(dot(xi(preN,:),xi(preN,:)));
    temxi(preN,:) = xi(preN,(setdiff((1:1:(Ny+1)*(Nx+1)),BC)));
    %% update AAL1 ALL2
    ALL1(preN,1:preN)=temxi(preN,:)*(CoeffMatrixL1*temxi(1:preN,:)');
    ALL1(1:preN,preN)=temxi(1:preN,:)*(CoeffMatrixL1*temxi(preN,:)');
    ALL2(preN,1:preN)=temxi(preN,:)*(CoeffMatrixL2*temxi(1:preN,:)');
    ALL2(1:preN,preN)=temxi(1:preN,:)*(CoeffMatrixL2*temxi(preN,:)');
end

if Nx==100
    save('al100.mat','temxi','xi','ALL1','ALL2')
elseif Nx==200
    save('al200.mat','temxi','xi','ALL1','ALL2')
elseif Nx==400
    save('al400.mat','temxi','xi','ALL1','ALL2')
else
    save('al800.mat','temxi','xi','ALL1','ALL2')
end

%% error tests
errL=zeros(1,Nmax);
lengthx=length(testeps);%length of eps
lengthy=length(testv);%length of potential v
length0=lengthx*lengthy;
urb=zeros(length0,(Ny+1)*(Nx+1));
w=zeros(TN,1);
%%
u_accurate=zeros(length0,(Ny+1)*(Nx+1));
k0=0;
for i1=1:lengthx
    for i2=1:lengthy
        k0=k0+1;
        v_test=testv(i2);
        eps_test=testeps(i1);
        u_accurate(k0,:)=FDMtwo(Nx, Ny,eps_test,v_test);
    end
end
BC=[(1:Nx+1:1+Ny*(Nx+1)),(Nx+1:Nx+1:(Ny+1)*(Nx+1))];
for iy=1:Ny+1
    F1(1,(iy-1)*(Nx-1)+1:(iy-1)*(Nx-1)+Nx-1)=gaussmf(xx,[0.1 0])*gaussmf(-1+hy*(iy-1),[0.1 0]);
end
%%
for N0=1:Nmax
    k0=0;
    for i1=1:lengthx
        eps_test=testeps(i1);
        lambdax=(eps_test/hx)^2;
        lambday =(eps_test/hy)^2;
        A1 = ALL1(1:N0,1:N0)*lambdax;
        A2 = ALL2(1:N0,1:N0)*lambday;
        for i2=1:lengthy
            k0=k0+1;
            v_test=testv(i2);
            for i3 = 1:Ny+1
                uold_test(1,(i3-1)*(Nx-1)+1:(i3-1)*(Nx-1)+Nx-1) = v_test*xx;
            end
            V_L = - v_test;
            V_R = v_test;
            errorPB=1;
            iter =1;
            while errorPB>10^(-8)
                BB = sparse((1:1:TN),(1:1:TN),...
                    cosh(uold_test),TN,TN);
                AD=temxi(1:N0,:)*(BB*temxi(1:N0,:)');
                %%
                B = A1+ A2 +AD;
                %% define right hand vector
                w(:,1)=cosh(uold_test).*uold_test-sinh(uold_test)+F1;
                w(1:(Nx-1):(Nx-1)*Ny+1,1)=w(1:(Nx-1):(Nx-1)*Ny+1,1)+lambdax*V_L;
                w(Nx-1:(Nx-1):TN)=w(Nx-1:(Nx-1):TN,1)+lambdax*V_R;
                F=temxi(1:N0,:)*w;
                %%
                coeff_test=B\F;
                % store all the nodes information of all paras
                temu=coeff_test'*temxi(1:N0,:);
                errorPB=max(abs((temu-uold_test)'));
                uold_test=temu;
                iter=iter+1;
            end
            urb(k0,setdiff((1:1:(Ny+1)*(Nx+1)),BC))=temu;
            urb(k0,1:Nx+1:1+Ny*(Nx+1))= V_L;% boundary condition
            urb(k0,Nx+1:Nx+1:(Ny+1)*(Nx+1))= V_R;% boundary condition
        end
    end
    ee=u_accurate-urb;
    errL(N0)=max(max(abs(ee)))/testv(lengthy);
    %         errL2(N0)=sqrt(max(sum((ee.^2),2)))/testv(lengthy);
    N0
end
u0=errL;
