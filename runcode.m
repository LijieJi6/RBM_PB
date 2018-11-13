
clear
tic
NN=[1000,2000,4000,8000];
% Nmax=18;% Select Nmax RB basis
% N=1000;
Paraeps=(0.08:0.02:0.4);
Parav=(0:0.25:5);
testeps=(0.085:0.01:0.395);
testv=(0.4:0.5:4.4);
%% offline and error test
% errL=zeros(4,Nmax);
% for in=1:4
%     N=NN(in);
% errL(in,:)=RBMone(Paraeps,Parav,N,Nmax,testeps,testv);
% end

%% online time 
% RBMonline


%% capacitance
N=10000;
Nmax=16;
epsc=0.08:0.005:0.4;
cap_exact=zeros(101,length(epsc));
for ic=1:length(epsc)
cap_exact(:,ic)=epsc(ic)*cosh((0:0.02:2)/2)/2;
end
[cap_rb,time_rb]=RBMcapacitance(Paraeps,Parav,N,Nmax,epsc);
Err=cap_rb-cap_exact;
err_abs=abs(Err);
%%
% mesh(epsc,0:0.02:2,err_abs,'LineWidth',3)
 mesh(epsc,0:0.02:2,Err,'LineWidth',3)
%  mesh(epsc,0:0.02:2,cap_rb,'LineWidth',3)
%  mesh(epsc,0:0.02:2,cap_exact,'LineWidth',3)
xlabel('\fontsize{15} D^{1/2}','FontName', 'Times New Roman')
ylabel('\fontsize{15} V','FontName', 'Times New Roman')
zlabel('\fontsize{18} Error','FontName', 'Times New Roman')
set(gca,'FontName','Times New Roman','FontSize',14,...
    'GridColor','k','FontWeight','bold','LineWidth',2)
box on
% colormap(cool)
%%
save('datacap.mat','cap_rb','cap_exact','err_abs','Err')
%%




