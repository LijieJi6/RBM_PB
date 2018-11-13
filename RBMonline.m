clear
 load offline
tic
preN=Nmax;
N=1000;
h=2/N;
B=zeros(N-1,preN);
w=zeros(N-1,1);
eps_test=0.1;
v_test=3.57;
lambda=eps_test^2/h^2;
AL0 = AL0*lambda;
xx=(-1+h:h:1-h);
u_old2=v_test'*xx;
temtemu(1)= -v_test;
temtemu(2)= +v_test;% boundary condition

error1=1;
while error1>10^(-8)
    BB =sparse((1:1:N-1),(1:1:N-1),...
        cosh(u_old2),N-1,N-1);
        AD=temxi*(BB*temxi');
    B = AL0 +AD;
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
toc
u_accurate=FDMPB(v_test,N,eps_test);
max(abs(u_accurate-u_test))
