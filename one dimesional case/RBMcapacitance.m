function [cap_RB,time_RB]=RBMcapacitance(Paraeps,Parav,N,Nmax,epsc)
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
save('alcap.mat','AL0','temxi','xi');

%% calculate capacitance
tic
h=2/N;
w=zeros(N-1,1);
xx=(-1+h:h:1-h);
Cap=zeros(101,length(epsc));
Q=zeros(101,1);
for ic=1:length(epsc)
eps_test=epsc(ic);
  lambda=eps_test^2/h^2;
    AL = AL0*lambda;
for i1= 0:101
    sigma_test= (i1-1)*0.02; 
    u_old2=sigma_test'*xx;
    temtemu(1)= -sigma_test;
    temtemu(2)= +sigma_test;% boundary condition
    error1=1;
    while error1>10^(-8)
        BB =sparse((1:1:N-1),(1:1:N-1),...
            cosh(u_old2),N-1,N-1);
        AD=temxi*(BB*temxi');
        B = AL +AD;
        %% define right hand vector
        w(1:N-1)=cosh(u_old2).*u_old2-sinh(u_old2);
        w(1)=w(1)+lambda*temtemu(1);
        w(N-1)=w(N-1)+lambda*temtemu(2);
        F=temxi*w;
        % solve this linear equation
        coeff_test=B\F;
        % store all the nodes information of all paras
        u_test=coeff_test'*temxi;
        error1=max(abs((u_test-u_old2)'));
        u_old2=u_test;
    end
    u_test=coeff_test'*xi(1:preN,:);
    Q(i1+1)=eps_test^2*(4*u_test(2)-3*u_test(1)-u_test(3))*N/4;
end
Cap(:,ic)=(Q(2:102)-Q(1:101))/0.02/2;% total capacitance
end
cap_RB=Cap;
time_RB=toc;




