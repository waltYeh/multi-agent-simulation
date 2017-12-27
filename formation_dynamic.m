function formation()
clear all
% L=[1,-1,0,0,0;
%     -1,3,-1,0,-1;
%     0,-1,3,-1,-1;
%     0,0,-1,2,-1;
%     0,-1,-1,-1,3];
% Delta=diag(diag(L));
% A=L-Delta;
t_final = 1;
t_steps = 11;
t = linspace(0,t_final,t_steps);
G=graph();
G=addnode(G,10);
ux=zeros(10,t_steps);
uy=zeros(10,t_steps);
ux(:,1)=[-1.2,-0.4,0.1,0.9,0.6]';
uy(:,1)=[-1.1,-0.4,0.6,0.9,0.1]';
G.Nodes.x=ux(:,1);
G.Nodes.y=uy(:,1);



R=1;
ksi_x=[R*sin(72*3/(180/pi)),R*sin(72*4/(180/pi)),0,R*sin(72/(180/pi)),R*sin(72*2/(180/pi))]';
ksi_y=[R*cos(72*3/(180/pi)),R*cos(72*4/(180/pi)),R,R*cos(72/(180/pi)),R*cos(72*2/(180/pi))]';
G.Nodes.ksi_x=ksi_x(:,1);
G.Nodes.ksi_y=ksi_y(:,1);
for i=1:length(t)
    T_next = t(i);
    clear Uode_x
    clear Uode_y
    if i~=1
        [Tode_x,Uode_x]=ode45(@(t,u)agreement_protocol_func(t,u,G,'x'),[t(i-1),T_next],ux(:,i-1));
        [Tode_y,Uode_y]=ode45(@(t,u)agreement_protocol_func(t,u,G,'y'),[t(i-1),T_next],uy(:,i-1));
        ux(:,i) = Uode_x(end,:);
        uy(:,i) = Uode_y(end,:);
    end
end
for i=1:11
%     figure(i)
    draw_graph(L,ux(:,i),uy(:,i),'b');
    pause(1)
    
end
end
function dx = agreement_protocol_func(t,x,G,axis)
dx = zeros(numnodes(G),1);
for i=1:numnodes(G)
    N=neighbors(G,i);
%     for j=1:N
        if axis=='x'
            dx(i)=dx(i)-(x(i)-x(j)-(G.Nodes.ksi_x(i)-G.Nodes.ksi_x(j)));
        elseif axis=='y'
            dx(i)=dx(i)-(x(i)-x(j)-(G.Nodes.ksi_y(i)-G.Nodes.ksi_y(j)));
        end
    end
end
end
function h=draw(L,ux,uy,arg)
hold on
n=length(L);
for i=1:n
   for j=1:n
       if i~=j
           if L(i,j)~=0
               h=plot([ux(i),ux(j)],[uy(i),uy(j)],arg);
           end
       end
   end
end
hold off
end