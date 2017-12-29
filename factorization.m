function factorization()
close all
s=[1 2 3];
t=[2 3 4];
G1=graph(s,t);
G2=graph(s,t);
% G3=graph(s,t);
% G1.Nodes.Position={[0,0,0],[1,0,0]}';
% G2.Nodes.Position={[0,0,0],[1,0,0]}';
% G3.Nodes.Position={[0,0,0],[1,0,0]}';
G=cartesian_product(G1,G2);
% G=cartesian_product(G12,G3);
figure(1)
plot(G)
L1=laplacian(G1);
L2=laplacian(G2);
L=laplacian(G);
t_final=4;

x1_0=[1;2;3;4];
y1_0=[1;1.2;1.4;1.6];
% z1_0=[0;0];
x2_0=[2;1.8;1.6;1.4];
y2_0=[1;2;3;4];
% z2_0=[0;0];
% x3_0=[2;3];
% y3_0=[3;2];
% z3_0=[0;1];

% Berechnung des dynamischen Ablaufs von einem einfachen Graph mit 2 Knoten
[Tode_x1,Uode_x1]=ode45(@(t,u)ode_func(t,u,L1),[0,t_final],x1_0);
[Tode_y1,Uode_y1]=ode45(@(t,u)ode_func(t,u,L1),[0,t_final],y1_0);
[Tode_x2,Uode_x2]=ode45(@(t,u)ode_func(t,u,L1),[0,t_final],x2_0);
[Tode_y2,Uode_y2]=ode45(@(t,u)ode_func(t,u,L1),[0,t_final],y2_0);
% [Tode_x3,Uode_x3]=ode45(@(t,u)ode_func(t,u,L1),[0,t_final],x3_0);
% [Tode_y3,Uode_y3]=ode45(@(t,u)ode_func(t,u,L1),[0,t_final],y3_0);
t_fix=0:0.2:t_final;
x1 = interp1(Tode_x1,Uode_x1,t_fix); 
y1 = interp1(Tode_y1,Uode_y1,t_fix); 
x2 = interp1(Tode_x2,Uode_x2,t_fix); 
y2 = interp1(Tode_y2,Uode_y2,t_fix); 
% x3 = interp1(Tode_x3,Uode_x3,t_fix); 
% y3 = interp1(Tode_y3,Uode_y3,t_fix); 
figure(3)
subplot(2,1,1)
plot(Tode_x1,Uode_x1)
title('Berechnung des dynamischen Ablaufs von einem einfachen Graph mit 2 Knoten')
title('Dynamics of atom graph G1')
xlabel('t(s)')
subplot(2,1,2)
plot(Tode_x2,Uode_x2)
title('Berechnung des dynamischen Ablaufs von einem einfachen Graph mit 2 Knoten')
title('Dynamics of atom graph G2')
xlabel('t(s)')
%Berechnung des Ablaufs von einem Kompositgraph
%aus Kronecker Product vom Ergebnis des einfachen Graphs
x=zeros(numnodes(G),length(t_fix));
y=zeros(numnodes(G),length(t_fix));
for i=1:length(t_fix)
    x(:,i)=kron(x1(i,:),x2(i,:));
    y(:,i)=kron(y1(i,:),y2(i,:));
end
% for i=1:length(Uode_y1)
%     y(:,i)=kron(kron(Uode_y1(i,:),Uode_y1(i,:)),Uode_y1(i,:));
% end
figure(4)
plot(t_fix,x)
title('Berechnung des Ablaufs von einem Kompositgraph aus Kronecker Product vom Ergebnis des einfachen Graphs')
title('Dynamics of composite graph of G1 and G2')
xlabel('t(s)')
% Berechnung des dynamischen Ablaufs von einem Kompositgraph aus
% Cartesian_product dreier einfachen Graphen
x_0=kron(x1_0,x2_0);
[Tode_x,Uode_x]=ode45(@(t,u)ode_func(t,u,L),[0,t_final],x_0);
y_0=kron(y1_0,y2_0);
[Tode_y,Uode_y]=ode45(@(t,u)ode_func(t,u,L),[0,t_final],y_0);
% [Tode_y1,Uode_y1]=ode45(@(t,u)ode_func(t,u,L1),[0,t_final],y1_0);
t_fix=0:0.2:t_final;
x = interp1(Tode_x,Uode_x,t_fix)'; 
y = interp1(Tode_y,Uode_y,t_fix)'; 
figure(5)
plot(Tode_x,Uode_x)
title('Berechnung des dynamischen Ablaufs von einem Kompositgraph aus Cartesian_product dreier einfachen Graphen')
figure(6)
draw_graph(L,x(:,1),y(:,1),'b-o')
hold on
draw_graph(L,x(:,9),y(:,9),'r-o')
draw_graph(L,x(:,14),y(:,14),'y-o')
hold off
axis equal

figure(2)
subplot(1,2,1)
draw(L1,x1_0,y1_0,'-o')
% hold on
% draw_graph(L,x1(9,:),y1(9,:),'r-o')
% draw_graph(L,x1(14,:),y1(14,:),'y-o')
% hold off
axis equal
axis([-1 5 0 5])
subplot(1,2,2)
draw(L2,x2_0,y2_0,'-o')
% hold on
% draw_graph(L,x2(9,:),y2(9,:),'r-o')
% draw_graph(L,x2(14,:),y2(14,:),'y-o')
% hold off
axis equal
axis([-1 5 0 5])
end
function G12=cartesian_product(G1,G2)
L1=laplacian(G1);
L2=laplacian(G2);
L12=kron(L1,eye(length(L2),length(L2)))+kron(eye(length(L1),length(L1)),L2);
A12=L12-diag(diag(L12));
G12=graph(A12);
end
function du = ode_func(t,u,L)
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