function decentralized_formation()
close all
global decent_err1
global decent_err2
global time
decent_err1=[];
decent_err2=[];
time=[];
goal=15+15i;
% obst=[];
obst=[3.5;2+3i;4+2i;3+3.5i;4+4.5i;5+3.2i;6+4.6i;5.2+6i;
    8+7.5i];
vdes=1;
k=[1;2;2;2];
% k_obst, k_v, k_x k_form
t_final = 30;
sd = [2 2 3];
td = [1 3 4];
zref1=[0.4+0.4i;-0.4-0.4i;-0.4-0.4i];%--
zref2=[0.6i-0.5;-1;0.5-0.6i];%z
EdgeTable = table([sd' td'],zref1,zref2,'VariableNames',{'EndNodes' 'zref1' 'zref2'});
D = digraph(EdgeTable);
r0=[0;-2-4.75i;0;0];
r0(1)=r0(2)+zref2(1);
r0(3)=r0(2)+zref2(2);
r0(4)=r0(3)+zref2(3);
v0=[0,0,0,0]';
u0=[v0;r0];
op=odeset('Events',@(t,y)eventfun(t,y,goal));
figure(1)
hold on
axis equal
plot(goal,'o')
plot(obst,'o')
grid
[Tode,Uode]=ode45(@(t,u)totalized_agent( t,u,D,k,goal,obst,vdes),[0,t_final],u0,op);
vel = Uode(:,1:numnodes(D));
pos = Uode(:,numnodes(D)+1:numnodes(D)*2);
err1=zeros(length(Tode),1);
err2=zeros(length(Tode),1);
for ii=1:length(Tode)
    [err1(ii),err2(ii)]=centralized_form_err(pos(ii,:),D);
end
delta_err=err1-err2;
zero_pts=pickzero(delta_err);
for ii=1:numnodes(D)
    plot(pos(:,ii))
end
for ii=1:length(zero_pts)
    plot(real(pos(zero_pts(ii),:)),imag(pos(zero_pts(ii),:)),'--*')
end
figure(2)
plot(Tode,abs(vel),'*')
figure(3)
subplot(2,1,1)
hold on
plot(time,decent_err1(:,2))
plot(time,decent_err1(:,1))
plot(time,decent_err1(:,3))
plot(time,decent_err1(:,4))

plot(time,decent_err2(:,2),'--')
plot(time,decent_err2(:,1),'--')
plot(time,decent_err2(:,3),'--')
plot(time,decent_err2(:,4),'--')

title('Global formation error consensus while flying towards target, avoiding obstacles, and maintaining either of 2 formations')
legend('form 1 agnt 1','form 1 agnt 2','form 1 agnt 3','form 1 agnt 4','form 2 agnt 1','form 2 agnt 2','form 2 agnt 3','form 2 agnt 4')
xlabel('time(s)')
ylabel('decentralized error estimation')
grid

subplot(2,1,2)
plot(Tode,err1,Tode,err2)
title('Global formation error while flying towards target, avoiding obstacles, and maintaining either of 2 formations')
xlabel('time(s)')
ylabel('centralized error estimation')
legend('form 1 (in line)','form 2 (Z form)')
grid
end


function du=totalized_agent(t,u,D,k,goal,obst,vmax)
global decent_err1
global decent_err2
global time
n=numnodes(D);
du=zeros(n*2,1);
v_nbrs=u(1:n);
x_nbrs=u(n+1:2*n);

loc_err1=zeros(n,1);
loc_err2=zeros(n,1);
for ii=1:n
    [loc_err1(ii),loc_err2(ii)]=local_form_err(x_nbrs,D,ii);
end
u0=[loc_err1;loc_err2];
[Tode,Uode]=ode45(@(t,u)totalized_err_syn(t,u,D),[0,3],u0);
glob_e1=Uode(end,1:n);
glob_e2=Uode(end,n+1:n*2);

decent_err1=[decent_err1;glob_e1];
decent_err2=[decent_err2;glob_e2];
time=[time;t];
for ii=1:n
    u_agent=zeros(2,1);
    u_agent(1)=u(ii);
    u_agent(2)=u(ii+n);
    [du_agent] = ...
        ctrl_law_agent( t,u_agent,D,k,goal,obst,vmax,ii,x_nbrs,v_nbrs,glob_e1,glob_e2 );
    du(ii)=du_agent(1);
    du(ii+n)=du_agent(2);
end
end

function du=totalized_err_syn(t,u,D)
n=numnodes(D);
e1=u(1:n);
e2=u(n+1:2*n);
du=zeros(2*n,1);
for ii=1:n
    [de1,de2] = err_syn_agent(D,ii,e1,e2 );
    du(ii)=de1;
    du(n+ii)=de2;
end
end

function [de1,de2] = err_syn_agent(D,nodeIndex,err1_nbrs,err2_nbrs )
N=predecessors(D,nodeIndex);
N=[N;successors(D,nodeIndex)];
de1=0;
de2=0;
for ii=1:length(N)
    de1=de1-(err1_nbrs(nodeIndex)-err1_nbrs(N(ii)));
end
for ii=1:length(N)
    de2=de2-(err2_nbrs(nodeIndex)-err2_nbrs(N(ii)));
end
end

function [du] = ctrl_law_agent( t,u,D,k,goal,obst,vmax,nodeIndex,x_nbrs,v_nbrs,err1_nbrs,err2_nbrs )
n=numnodes(D);
du=zeros(2,1);
n_obst=length(obst);
pos=u(2);
vel=u(1);
if(err1_nbrs(nodeIndex)>err2_nbrs(nodeIndex))
    zref=D.Edges.zref2;
    agreement=2;
else
    zref=D.Edges.zref1;
    agreement=1;
end
form_err=0;
d_form_err=0;
for jj=1:numnodes(D)
    EdgeIdx=findedge(D,jj,nodeIndex);
    if EdgeIdx~=0
        form_err=form_err+zref(EdgeIdx)-(pos-x_nbrs(jj));
        d_form_err=d_form_err+(vel-v_nbrs(jj));
    end
end
du(1)=du(1)+k(4)*form_err;%-k(4)*d_form_err;
for jj=1:n_obst
    du(1)=du(1)-k(1)*(obst(jj)-pos)/abs(obst(jj)-pos)^3;
end
x_err=goal-pos;
vdes=k(3)*x_err;
norm_vdes=abs(vdes);
if norm_vdes>vmax
    vdes=vmax*vdes/norm_vdes;
end
v_err=vdes-vel;
du(1)=du(1)+k(2)*v_err;
du(2)=vel;
end


function [err1,err2]=centralized_form_err(x,D)
% n=numnodes(D);
err1=0;
err2=0;
for ii=1:numedges(D)
    err1=err1+(abs( x(D.Edges.EndNodes(ii,1)) - x(D.Edges.EndNodes(ii,2)) )^2-( abs(D.Edges.zref1(ii)) )^2)^2;
    err2=err2+(abs( x(D.Edges.EndNodes(ii,1)) - x(D.Edges.EndNodes(ii,2)) )^2-( abs(D.Edges.zref2(ii)) )^2)^2;
end
end
function [err1,err2]=local_form_err(x,D,nodeIndex)
err1=0;
err2=0;
for jj=1:numnodes(D)
    EdgeIdx=findedge(D,nodeIndex,jj);
    if EdgeIdx==0
        EdgeIdx=findedge(D,jj,nodeIndex);
    end
    if EdgeIdx~=0
        err1=err1+(abs( x(D.Edges.EndNodes(EdgeIdx,1)) - x(D.Edges.EndNodes(EdgeIdx,2)) )^2-( abs(D.Edges.zref1(EdgeIdx)) )^2)^2;
        err2=err2+(abs( x(D.Edges.EndNodes(EdgeIdx,1)) - x(D.Edges.EndNodes(EdgeIdx,2)) )^2-( abs(D.Edges.zref2(EdgeIdx)) )^2)^2;
    end
end
end

function [value,isterminal,direction] = eventfun(t,y,goal)
n2=length(y);
n=n2/2;
avr_pos=sum(y(n+1:n2))/n;
dist=abs(avr_pos-goal);
value = dist-0.1;
isterminal= 1;
direction = 0;
end

function zeroindex=pickzero(x) 
%pick the index of 0 or the nearest in a vector  
m = length(x); 
x1=x(1:m-1); 
x2=x(2:m); 
indz = find(x==0); %zero point 
indzer = find(x1.*x2<0); %negative/positive 
n=length(indzer); 
for i=1:n 
	if abs(x(indzer(i)))>abs(x(indzer(i)+1)) 
        indzer(i)=indzer(i)+1; 
	end 
end 
zeroindex=sort([indz indzer],2,'ascend');
end