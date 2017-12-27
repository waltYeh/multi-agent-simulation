function centralized_formation()
goal=12+12i;
% obst=[];
obst=[3+2i;2+3i;4+3i;3+4i;5+5i];
vdes=1;
k=[1;2;1;5];
% k_obst, k_v, k_x k_form
t_final = 30;
s = [1 2 3];
t = [2 3 1];
G = graph(s,t);
sd = [2 2];
td = [1 3];
D = digraph(sd,td);
incidence_D=[1,0;-1,-1;0,1];
% zref1=[1+1i;-1-1i]/2;
zref1=[0.5+0.5*sqrt(3)*1i;0.5-sqrt(3)/2+sqrt(3)/2-0.5i];

r0=[1+1i;1;0];
v0=[0.2i,0.2i,0.2i]';
u0=[v0;r0];
[Tode,Uode]=ode45(@(t,u)ctrl_law( t,u,D,k,goal,obst,vdes,zref1,[0;0] ),[0,t_final],u0);
vel = Uode(:,1:3);
pos = Uode(:,4:6);

figure(1)
hold on
for i=1:3
    plot(pos(:,i))
%     quiver(real(pos(:,i)),imag(pos(:,i)),cos(theta(:,i)),sin(theta(:,i)))
%     quiver(real(pos(1:10:end,i)),imag(pos(1:10:end,i)),real(vel(1:10:end,i)),imag(vel(1:10:end,i)))
    axis equal
end
for i=1:20:length(Tode)
    plot(real(pos(i,:)),imag(pos(i,:)),'--')
end
plot(goal,'+')
plot(obst,'+')
hold off
grid
figure(2)
plot(Tode,abs(vel))
% for i=1:10
%     plot(Tode,real(Uode(:,i+10)))
%     hold on
% end
% hold off
% grid
% figure(3)
% for i=1:10
%     plot(Tode,imag(Uode(:,i+10)))
%     hold on
% end
% hold off
% grid
% figure(6)
% for i=1:10
% plot(r(i,:));
% hold on
% end
% hold on
% plot(CoM,'-*');
% draw_graph(laplacian(G),real(r(:,end)),imag(r(:,end)),'-+b')
% hold off
% axis equal
% grid
end



function du = ctrl_law( t,u,D,k,goal,obst,vmax,zref,d_zref )
n=numnodes(D);
G=graph(D.Edges);
LD=laplacian(G);
% invD=flipedge(D);
% LD_inv=diag(indegree(invD))-adjacency(invD);
% LD=diag(indegree(D))-adjacency(D);
du=zeros(n*2,1);
n_obst=length(obst);
x=u(n+1:2*n);
vel=u(1:n);
a_form=-k(4)*LD*x-k(4)*LD*vel+k(4)*incidence(D)*zref+k(4)*incidence(D)*d_zref;
for ii=1:n
%     F_obst
    for jj=1:n_obst
        du(ii)=du(ii)-k(1)*(obst(jj)-x(ii))/abs(obst(jj)-x(ii))^3;
    end
%     F_goal
    x_err=goal-x(ii);
    vdes=k(3)*x_err;
    norm_vdes=abs(vdes);
    if norm_vdes>vmax
        vdes=vmax*vdes/norm_vdes;
    end
    v_err=vdes-vel(ii);
    du(ii)=du(ii)+k(2)*v_err;
end
du(1:n)=du(1:n)+a_form;
for ii=1:n
    du(ii+n)=vel(ii);
end
end
