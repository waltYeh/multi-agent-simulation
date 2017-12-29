function centralized_formation()
close all
goal=12+12i;
% obst=[];
obst=[3+2i;2+3i;4+3i;3+4i;5+4.8i;5+6.5i];
vdes=1;
k=[1;2;2;2];
% k_obst, k_v, k_x k_form
t_final = 30;
sd = [2 2 3];
td = [1 3 4];
zref1=[1+1i;-1-1i;-1-1i]/2;
zref2=[1i;-1;-1i]/2;
EdgeTable = table([sd' td'],zref1,zref2,'VariableNames',{'EndNodes' 'zref1' 'zref2'});
D = digraph(EdgeTable);
% incidence_D=[1,0;-1,-1;0,1];



r0=[1+1i;0.5+0.5i;0;1];
r0=r0-3-3i;
v0=[0,0,0,0]';
u0=[v0;r0];
op=odeset('Events',@(t,y)eventfun(t,y,goal));
[Tode,Uode]=ode45(@(t,u)ctrl_law( t,u,D,k,goal,obst,vdes),[0,t_final],u0,op);
vel = Uode(:,1:numnodes(D));
pos = Uode(:,numnodes(D)+1:numnodes(D)*2);
err1=zeros(length(Tode),1);
err2=zeros(length(Tode),1);
for ii=1:length(Tode)
    [err1(ii),err2(ii)]=centralized_form_err(pos(ii,:),D);
end
figure(1)
hold on
axis equal
plot(goal,'o')
plot(obst,'o')
grid
for i=1:20:length(Tode)
    plot(real(pos(i,:)),imag(pos(i,:)),'--')
end
for i=1:numnodes(D)
    plot(pos(:,i))
%     quiver(real(pos(:,i)),imag(pos(:,i)),cos(theta(:,i)),sin(theta(:,i)))
%     quiver(real(pos(1:10:end,i)),imag(pos(1:10:end,i)),real(vel(1:10:end,i)),imag(vel(1:10:end,i)))
end
hold off
figure(2)
plot(Tode,abs(vel),'*')
figure(3)
plot(Tode,err1,Tode,err2)
end



function du = ctrl_law( t,u,D,k,goal,obst,vmax )
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
[err1,err2]=centralized_form_err(x,D);
if(err1>err2)
    zref=D.Edges.zref2;
else
    zref=D.Edges.zref1;
end
d_zref=zeros(numedges(D),1);
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
function [err1,err2]=centralized_form_err(x,D)
% n=numnodes(D);
% G=graph(D.Edges);
err1=0;
err2=0;
for ii=1:numedges(D)
    err1=err1+(abs( x(D.Edges.EndNodes(ii,1)) - x(D.Edges.EndNodes(ii,2)) )^2-( abs(D.Edges.zref1(ii)) )^2)^2;
    err2=err2+(abs( x(D.Edges.EndNodes(ii,1)) - x(D.Edges.EndNodes(ii,2)) )^2-( abs(D.Edges.zref2(ii)) )^2)^2;
end
end
function [value,isterminal,direction] = eventfun(t,y,goal)
n2=length(y);
n=n2/2;
avr_pos=sum(y(n+1:n2))/n;
dist=abs(avr_pos-goal);
value = dist-1;
isterminal= 1;
direction = 0;
end