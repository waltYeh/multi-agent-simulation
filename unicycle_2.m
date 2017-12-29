function unicycle_2()
m=1;
rate0=0.2;
k=1;
t_final = 30;
% s = [1 1 1 2 2 3 3 4 5 5 6 7 7 8 1  2 9];
% t = [2 4 8 3 7 4 6 5 6 8 7 8 9 9 10 10 10];
s = [1 2 3 4 5 6 7 8 9 10];
t = [2 3 4 5 6 7 8 9 10 1];
G = graph(s,t);
% G=graph();
% G=addnode(G,10);
r0=(0:9)'+1i*(9:-1:0)';
theta0=[1,2,-1,0,-1.5,1.2,2.2,-1.2,0.2,-1.7]';
% theta0=[1,2,-1,0,-1.5,1,2,-1,0,-1.5]'/10;

% if theta0 are identical, no balancing will occure when no disturbance of
% yaw angle comes to present
% k positive, balanced; k negative, synchronized
[Tode,Uode]=ode45(@(t,u)part_connect_ctrl_fun( t,u,G,k,m,rate0 ),[0,t_final],[theta0,r0]);
theta = Uode(:,1:10)';
r = Uode(:,11:20)';
states=zeros(1,length(Tode));
potentials=zeros(1,length(Tode));
CoM=zeros(1,length(Tode));
for i=1:length(Tode)
    states(i)=state(m,theta(:,i),G);
    potentials(i)=potential(m,theta(:,i),G);
    CoM(i)=sum(r(:,i))/10;
end
figure(1)
for i=1:10
    plot(Tode,Uode(:,i))
    hold on
end
hold off
grid
figure(2)
for i=1:10
    plot(Tode,real(Uode(:,i+10)))
    hold on
end
hold off
grid
figure(3)
for i=1:10
    plot(Tode,imag(Uode(:,i+10)))
    hold on
end
hold off
grid
figure(4)
plot(Tode,real(states),Tode,imag(states))
grid
figure(5)
plot(Tode,real(potentials),Tode,imag(potentials))
grid
figure(6)
for i=1:10
plot(r(i,:));
hold on
end
hold on
%     h2=plot(CoM,'-*');

% for i=1:length(r)
%     h1(i)=plot(r(:,i)','o');
%     h2(i)=plot(CoM(i)','*');
%     axis([-50,50,-50,50])
% %     plot(CoM(i)','--*')
%     hold on
%     if i-4>0
% 
%         set(h1(i-4),'visible','off');
%         set(h2(i-4),'visible','off');
%     end
%     pause(0.1)
% end
% plot(r(:,end),'-+')
draw_graph(laplacian(G),real(r(:,end)),imag(r(:,end)),'-+b')
hold off
axis equal
grid
end

function res=state(m,theta,G)
n=numnodes(G);
res=1/n/m*ones(n,1)'*exp(1i*m*theta);
end

function res=potential(m,theta,G)
n=numnodes(G);
res=1/2/n/m^2*conj(exp(1i*m*theta))'*ones(n,1)*ones(n,1)'*exp(i*m*theta);
end

function res=critical_pts(m,theta,G)
res=1/2*conj(exp(1i*m*theta))'*laplacian(G)*exp(i*m*theta);
end

function res=q(rate0,theta,G,r)
n=numnodes(G);
res=zeros(n,1);
for i=1:n
    res(i)=exp(1i*theta(i))-1i*rate0*r(i);
end
end

function du = ode_ctrl_fun( t,u,G,k,m,rate0 )
n=numnodes(G);
du=zeros(n*2,1);
for ii=1:n
    for j=1:n
        du(ii)=du(ii)-sin(u(j)-u(ii));
    end
    du(ii)=du(ii)*k/n+rate0;
end
for ii=1+n:n+n
    du(ii)=exp(1i*u(ii-10));
end
end

function du = part_connect_ctrl_fun( t,u,G,k,m,rate0 )
n=numnodes(G);
du=zeros(n*2,1);
for ii=1:n
    N=neighbors(G,ii);
    for j=1:length(N)
        du(ii)=du(ii)-sin(m*(u(N(j))-u(ii)));
    end
    du(ii)=du(ii)*k/n*m+rate0;
end
for ii=1+n:n+n
    du(ii)=exp(1i*u(ii-10));
end
end

function du = grad_flow_ctrl_fun( t,u,G,k,m,rate0 )
n=numnodes(G);
du=zeros(n*2,1);
L=laplacian(G);
qs=q(rate0,u(1:10),G,u(11:20));
for ii=1:n
    a=L(ii,:)*qs;
    b=1i*exp(1i*u(ii));
%     du(ii)=k*dot(a,1i*exp(1i*u(ii)))+rate0;
    du(ii)=k*(real(a)*real(b)+imag(a)*imag(b))+rate0;
end
for ii=1+n:n+n
    du(ii)=exp(1i*u(ii-10));
end
end