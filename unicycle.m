function unicycle()
t_final = 10;
t_steps = 11;
t = linspace(0,t_final,t_steps);
G=graph();
G=addnode(G,10);
r=zeros(10,t_steps);
r(:,1)=(0:9)'+1i*(9:-1:0)';
theta=zeros(10,t_steps);
states=zeros(1,t_steps);
potentials=zeros(1,t_steps);
theta(:,1)=[1,2,-1,0,-1.5,1.2,2.2,-1.2,0.2,-1.7]';
% theta(:,1)=ones(10,1)'*5;

T_all=zeros(1,1);
U_all=zeros(10,1);
% u=zeros(10,t_steps);
% sta=state(1,G);
% pot=potential(1,G);
for i=1:length(t)
    T_next = t(i);
    clear Uode
    clear Tode
%     clear Uode_y
    if i~=1
        [Tode,Uode]=ode45(@(t,u)ode_ctrl_fun( t,u,G,-1 ),[t(i-1),T_next],theta(:,i-1));
%         [Tode_y,Uode_y]=ode45(@(t,u)ode_ctrl_fun( t,u,G,k ),[t(i-1),T_next],uy(:,i-1));
        theta(:,i) = Uode(end,:);
%         uy(:,i) = Uode_y(end,:);
        T_all=[T_all,Tode'];
        U_all=[U_all,Uode'];
    end
end
for i=1:length(t)
    states(i)=state(1,theta(:,i),G);
    potentials(i)=potential(1,theta(:,i),G);
end
for i=1:10
    plot(T_all,U_all(i,:))
    hold on
end
hold off
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
du=zeros(n,1);
for ii=1:n
    for j=1:n
        du(ii)=du(ii)-sin(u(j)-u(ii));
    end
    du(ii)=du(ii)*k/n;
%     du(ii,2)=exp(1i*du(ii,1));
end
end