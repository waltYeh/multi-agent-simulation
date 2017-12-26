function unicycle_2()
t_final = 10;
G=graph();
G=addnode(G,10);
r0=(0:9)'+1i*(9:-1:0)';
theta0=[1,2,-1,0,-1.5,1.2,2.2,-1.2,0.2,-1.7]'/10;
theta0=[1,2,-1,0,-1.5,1,2,-1,0,-1.5]'/10;

% if theta0 are identical, no balancing will occure when no disturbance of
% yaw angle comes to present
% k positive, balanced; k negative, synchronized
[Tode,Uode]=ode45(@(t,u)ode_ctrl_fun( t,u,G,1 ),[0,t_final],[theta0,r0]);
theta = Uode(:,1:10)';
r = Uode(:,11:20)';
states=zeros(1,length(Tode));
potentials=zeros(1,length(Tode));
CoM=zeros(1,length(Tode));
for i=1:length(Tode)
    states(i)=state(1,theta(:,i),G);
    potentials(i)=potential(1,theta(:,i),G);
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
plot(r')
hold on
plot(CoM','--*')
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
function du = ode_ctrl_fun( t,u,G,k )
n=numnodes(G);
du=zeros(n*2,1);
for ii=1:n
    for j=1:n
        du(ii)=du(ii)-sin(u(j)-u(ii));
    end
    du(ii)=du(ii)*k/n;
%     du(ii,2)=exp(1i*du(ii,1));
end
for ii=1+n:n+n
    du(ii)=exp(1i*u(ii-10));
end
    
end