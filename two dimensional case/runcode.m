
%% This is the main code
clear
tic
NN=[100,200,400,800];
Nmax=20;% Select Nmax RB basis
Paraeps=(0.08:0.02:0.4);
Parav=(0:0.25:5);
testeps=(0.085:0.01:0.395);
testv=(0.4:0.5:4.4);
%% offline and error test
% errL=zeros(4,Nmax);
% for in=1:4
%     Nx=NN(in);
%     errL(in,:)=RBMtwo(Paraeps,Parav,Nx,Nmax,testeps,testv);
% end


%% calculate differential capacitance by FDM or RBM
Nx=400;
% load al200 % al200 saves the needed offline data
epsc=0.08:0.005:0.4;
% epsc=0.08:0.01:0.4;
% [cap_fd,timefdm]=cap_fdm(Nx,epsc);% capacitance by fdm
% [cap_rb,timerbm]=cap_rbm(Nx,epsc);% capacitance by rbm
Err=cap_rb-cap_fd;
err_abs=abs(Err);


% save ('cap.mat','cap_fd','cap_rb','Err','timerbm','timefdm')
% %%
% % mesh(epsc,0:0.02:2,err_abs,'LineWidth',3)
%  mesh(epsc,0:0.02:2,Err,'LineWidth',3)
 mesh(epsc,0:0.02:2,cap_rb,'LineWidth',3)
% %  mesh(epsc,0:0.02:2,cap_exact,'LineWidth',3)
xlabel('\fontsize{15} D^{1/2}','FontName', 'Times New Roman')
ylabel('\fontsize{15} V','FontName', 'Times New Roman')
zlabel('\fontsize{15} Differential Capacitance','FontName', 'Times New Roman')
set(gca,'FontName','Times New Roman','FontSize',14,...
    'GridColor','k','FontWeight','bold','LineWidth',2)
box on
% % colormap(cool)
% %%
% save('datacap.mat','cap_rb','cap_exact','err_abs','Err')
% %%




