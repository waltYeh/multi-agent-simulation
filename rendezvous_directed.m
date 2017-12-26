function rendezvous_directed()
Adir=[0,1,0,0,0;
    0,0,0,0,1;
    0,1,0,0,1;
    0,0,1,0,0;
    0,0,0,-1,0];
Deltadir=diag(Adir*ones(5,1));
L=Deltadir-Adir;
t_final = 4;
t_steps = 41;
t = linspace(0,t_final,t_steps);
ux=zeros(t_steps,4);
uy=zeros(t_steps,4);
ux(1,:)=[0,-1,2,3];
uy(1,:)=[2,1,0,2];
figure
hold on
draw(L,ux(1,:),uy(1,:),'k')
for i=2:length(t)
    T_next = t(i);
    clear Uode_x
    clear Uode_y
    [Tode_x,Uode_x]=ode45(@(t,u)ode_func(t,u,L),[0,T_next],ux(1,:));
    [Tode_y,Uode_y]=ode45(@(t,u)ode_func(t,u,L),[0,T_next],uy(1,:));
    ux(i,:) = Uode_x(end,:);
    uy(i,:) = Uode_y(end,:);
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