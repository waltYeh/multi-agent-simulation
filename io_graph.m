function io_graph()
%network partition
% Df=[-1,0,0,1,0,0,-1,0;
%     0,0,1,-1,-1,0,0,0;
%     1,-1,0,0,0,0,0,1;
%     0,1,-1,0,0,1,0,0];
% Di=[0,0,0,0,1,-1,0,0;
%     0,0,0,0,0,0,1,-1];
Df=[-1,0,0;
    1,-1,0;
    0,1,-1];
Di=[0,0,1];
D=[Df;Di];
Af=Df*Df';
Ai=Di*Di';
Bf=Df*Di';
L=[Af,Bf;Bf',Ai];
% L=D*D';
Adj=L-diag(diag(L));
G=graph(Adj);
[eigvec,eigval] = eig(Af);
% for i=1:size(eigvec,2)
%     for j=1:size(Bf,2)
%         dot(eigvec(:,i),Bf(:,j))
%     end
% end
% controllable and observable if and only if
% none of eigvector(Af) are simultaneously 
% orthogonal to all columns of Bf
t_final = 4;
xf0=[0,1,2];
yf0=[0,2,0];
sys=ss(-Af,-Bf,-Bf',0);
[u,t] = gensig('square',4,6,0.1);
[Uode_yx,Uode_tx,Uode_xf]=lsim(sys,u,t,xf0);
[Uode_yy,Uode_ty,Uode_yf]=lsim(sys,u,t,yf0);
plot(Uode_xf,Uode_yf)
% [Tode_x,Uode_x]=ode45(@(t,u)ode_func(t,u,Af),[0,t_final],ux0);
% plot(G)
end
function du = ode_func(t,u,L)
% L = [-2,1,1,0;
%     1,-3,1,1;
%     1,1,-3,1;
%     0,1,1,-2];
du = -L*u;
end