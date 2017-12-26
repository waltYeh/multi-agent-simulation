function rendezvous_undirected()
AdjHalf=[0,1,1,0;
    0,0,1,0;
    0,0,0,1;
    0,0,0,0];
Adj=AdjHalf'+AdjHalf;
G=graph(Adj);
figure(1)
plot(G)
figure(2)
Deg=zeros(length(Adj),length(Adj));
for i=1:length(Adj)
    Deg(i,i)=sum(Adj(:,i));
end
L=Deg-Adj;
Incidence=[1,-1,0,0;
    0,1,-1,0;
    -1,0,1,-1;
    0,0,0,1];
Le=Incidence'*Incidence;
null(Le)
null(L)
t_final = 4;
t_steps = 41;
t = linspace(0,t_final,t_steps);
ux=zeros(t_steps,4);
uy=zeros(t_steps,4);
ex=zeros(t_steps,4);
ey=zeros(t_steps,4);
ex_check=zeros(t_steps,4);
ey_check=zeros(t_steps,4);
edge_disagree_x=zeros(t_steps,1);
edge_disagree_y=zeros(t_steps,1);
vertex_disagree_x=zeros(t_steps,1);
vertex_disagree_y=zeros(t_steps,1);
check_disagr_ineq=zeros(t_steps,2);
ux(1,:)=[0,-1,2,3];
uy(1,:)=[2,1,0,2];
ex(1,:)=Incidence'*ux(1,:)';
ey(1,:)=Incidence'*uy(1,:)';
ex_check(1,:)=Incidence'*ux(1,:)';
ey_check(1,:)=Incidence'*uy(1,:)';
figure
hold on
draw(L,ux(1,:),uy(1,:),'k')
for i=1:length(t)
    T_next = t(i);
    clear Uode_x
    clear Uode_y
    if i~=1
        [Tode_x,Uode_x]=ode45(@(t,u)ode_func(t,u,L),[0,T_next],ux(1,:));
        [Tode_y,Uode_y]=ode45(@(t,u)ode_func(t,u,L),[0,T_next],uy(1,:));
        [Tode_x,Eode_x]=ode45(@(t,u)ode_func(t,u,Le),[0,T_next],ex(1,:));
        [Tode_y,Eode_y]=ode45(@(t,u)ode_func(t,u,Le),[0,T_next],ey(1,:));
        ux(i,:) = Uode_x(end,:);
        uy(i,:) = Uode_y(end,:);
        ex(i,:) = Eode_x(end,:);
        ey(i,:) = Eode_y(end,:);
    end
    ex_check(i,:)=Incidence'*ux(i,:)';
    ey_check(i,:)=Incidence'*uy(i,:)';
    vertex_disagree_x(i)=norm(ux(i,:));
    vertex_disagree_y(i)=norm(uy(i,:));
    edge_disagree_x(i)=norm(ex(i,:));
    edge_disagree_y(i)=norm(ey(i,:));
    check_disagr_ineq(i,1)=(edge_disagree_x(i)<=norm(Incidence)*(vertex_disagree_x(i)));
    check_disagr_ineq(i,2)=(edge_disagree_y(i)<=norm(Incidence)*(vertex_disagree_y(i)));
end
for i=[1,2,4,8,12,16,40]
    draw(L,ux(i,:),uy(i,:),'b')
end
end
function du = ode_func(t,u,L)
% L = [-2,1,1,0;
%     1,-3,1,1;
%     1,1,-3,1;
%     0,1,1,-2];
du = -L*u;
end
function draw(L,ux,uy,arg)
hold on
n=length(L);
for i=1:n
   for j=1:n
       if i~=j
           if L(i,j)~=0
               plot([ux(i),ux(j)],[uy(i),uy(j)],arg)
           end
       end
   end
end

end