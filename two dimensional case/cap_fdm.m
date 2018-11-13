function [u0,timefdm]=cap_fdm(Nx,epsc)
tic
Ny=Nx;
hx=2/Nx;
xx=-1+hx:hx:1-hx;
hy=hx;
for i3 = 1:Ny+1
    F1(1,(i3-1)*(Nx-1)+1:(i3-1)*(Nx-1)+Nx-1)=gaussmf(xx,[0.1 0])*gaussmf(-1+hy*(i3-1),[0.1 0]);
end
qq_fdm=zeros(102,Ny);
Cap_fdm=zeros(101,length(epsc));
for ic=1:length(epsc)
    eps_test=epsc(ic);
    for i1= 0:101
        v_test= (i1-1)*0.02;
        qq_fdm(i1+1,:)=capFDM(Nx,eps_test,v_test,F1);
    end
    Q_fdm=sum(qq_fdm,2);
Cap_fdm(:,ic)=(Q_fdm(2:102)-Q_fdm(1:101))/0.02/2;
end
timefdm=toc;
u0=Cap_fdm;


