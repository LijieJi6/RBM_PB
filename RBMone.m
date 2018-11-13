function u0=RBMone(Paraeps,Parav,N,Nmax,testeps,testv)
%% offline process and error tests

%Define variables:
%Paraeps -- train set
%Preeps(i) -- i th selected eps
%preN -- number of selected parameter
%xx -- variable x in (-1,1)
%h -- space step
%N -- partition number
h=2/N;
xx=(-1:2/N:1);
lengthx=length(Paraeps);%length of eps
lengthy=length(Parav);%length of potential v
preeps=zeros(1,Nmax);
prev=zeros(1,Nmax);

preeps(1)=0.24;% parameter center
prev(1)=2.5;
length0=lengthx*lengthy;
%% The first basis
xi(1,:)=FDMPB(prev(1),N,preeps(1));
xi(1,:)=xi(1,:)/sqrt(dot(xi(1,:),xi(1,:)));
temxi(1,:)=xi(1,2:N);
%%
eigval=zeros(1,length0);% store eigenvalues
errRB=zeros(1,Nmax);% error estimator of the selection
temu=zeros(length0,N+1);
temtemu=zeros(1,2);% boundary conditions
w=zeros(N-1,1);
mf=zeros(length0,N-1);
%%
CoeffMatrixL = sparse([1:1:N-1,1:1:N-2 ,2:1:N-1],[1:1:N-1, 2:1:N-1, 1:1:N-2],...
    [ 2*ones(1,N-1),-ones(1,N-2) ,-ones(1,N-2)],N-1,N-1);
ALL(1,1)=temxi(1,:)*(CoeffMatrixL*temxi(1,:)');
eig0 = eigs(CoeffMatrixL'*CoeffMatrixL,1,'sm');
for i1=1:lengthx
    % initial values for iteration
    u_old((i1-1)*lengthy+1:(i1-1)*lengthy+lengthy,:)=Parav'*xx;
    eigval((i1-1)*lengthy+1:(i1-1)*lengthy+lengthy)=Paraeps(i1)^4/h^4*eig0;
end
preN=1;
while preN<Nmax
    k0=0;
    for ix=1:lengthx
        lambda=Paraeps(ix)^2/h^2;
        for iy=1:lengthy
            iter=0;
            errPB=1;
            k0=k0+1;
            temtemu(1)= -Parav(iy);
            temtemu(2)= +Parav(iy);% boundary condition
            while errPB>10^(-8)
                AL(1:preN,1:preN)=lambda*ALL(1:preN,1:preN);
                % update diagonal
                BB =sparse((1:1:N-1),(1:1:N-1),cosh(u_old(k0,2:N)),N-1,N-1);
                AD(1:preN,1:preN)=temxi(1:preN,:)*(BB*temxi(1:preN,:)');
                B = AL(1:preN,1:preN) +AD(1:preN,1:preN);
                % right hand vector
                w(1:N-1)=cosh(u_old(k0,2:N)).*u_old(k0,2:N)-sinh(u_old(k0,2:N));
                w(1)=w(1)+lambda*temtemu(1);
                w(N-1)=w(N-1)+lambda*temtemu(2);
                F=temxi(1:preN,:)*w;
                % solve this linear equation
                coeff=B\F;
                % store all information of all paras
                temu(k0,:)=coeff'*xi(1:preN,:);
                errPB=max(abs(temu(k0,:)-u_old(k0,:)));
                u_old(k0,:)=temu(k0,:);
                iter=iter+1;
            end
            EB =lambda*CoeffMatrixL+BB;
            mf(k0,:) = EB*temu(k0,2:N)'-w;% residual vector of each parameter
        end
    end
    delta=sqrt(sum((mf.^2),2))./sqrt(eigval)';% L2 norm£¬£¬£¬
    [x,y]=max(delta);
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
    %% select a new basis
    xi(preN,:)=FDMPB(Parav(i2),N,Paraeps(i1));
    % G-S
    sum0=0;
    for m2=1:preN-1
        sum0=sum0+dot(xi(preN,:),xi(m2,:))/dot(xi(m2,:),xi(m2,:))*xi(m2,:);
    end
    xi(preN,:)=xi(preN,:)-sum0;
    xi(preN,:)=xi(preN,:)/sqrt(dot(xi(preN,:),xi(preN,:)));% normalization
    temxi(preN,:)=xi(preN,2:N);
    %% update matrix ALL
    ALL(preN,1:preN)=temxi(preN,:)*(CoeffMatrixL*temxi(1:preN,:)');
    ALL(1:preN,preN)=temxi(1:preN,:)*(CoeffMatrixL*temxi(preN,:)');
end
AL0 =ALL(1:preN,1:preN);
% save('offline.mat','AL0','temxi','xi','preeps','prev','preN');

%% test error
errL=zeros(1,Nmax);
lengthx=length(testeps);%length of eps
lengthy=length(testv);%length of potential v
length0=lengthx*lengthy;
u_accurate=zeros(length0,N+1);
urb=zeros(length0,N+1);
w=zeros(N-1,1);
%% high-fidelity solutions
k0=0;
for i1=1:lengthx
    for i2=1:lengthy
        k0=k0+1;
        v_test=testv(i2);
        eps_test=testeps(i1);
        u_accurate(k0,:)=FDMPB(v_test,N,eps_test);
    end
end

%% error tests
for N0=1:Nmax
    k0=0;
    for i1=1:lengthx
        eps_test=testeps(i1);
        lambda=eps_test^2/h^2;
        for i2=1:lengthy
            k0=k0+1;
            v_test=testv(i2);
            temtemu(1)= -v_test;
            temtemu(2)= +v_test;% boundary condition
            u_old2(1,:)=v_test'*xx;
            error1=1;
            while error1>10^(-8)
                BB =sparse((1:1:N-1),(1:1:N-1),...
                    cosh(u_old2(2:N)),N-1,N-1);
                AD=temxi(1:N0,:)*(BB*temxi(1:N0,:)');
                B = AL0(1:N0,1:N0)*lambda +AD;
                %% define right hand vector
                w(1:N-1)=cosh(u_old2(2:N)).*u_old2(2:N)-sinh(u_old2(2:N));
                w(1)=w(1)+lambda*temtemu(1);
                w(N-1)=w(N-1)+lambda*temtemu(2);
                F=temxi(1:N0,:)*w;
                % solve this linear equation
                coeff_test=B\F;
                % store all the nodes information of all paras
                u_test=coeff_test'*xi(1:N0,:);
                error1=max(abs((u_test-u_old2)'));
                u_old2=u_test;
            end
            urb(k0,:)=u_test;
        end
    end
    N0
    ee=u_accurate-urb;
    errL(N0)=max(max(abs(ee)))/testv(lengthy);% relative error
end
u0=errL;










