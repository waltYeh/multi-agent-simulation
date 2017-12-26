function factorization()
s=[1];
t=[2];
G1=graph(s,t);
G2=graph(s,t);
G3=graph(s,t);
% G1.Nodes.Position={[0,0,0],[1,0,0]}';
% G2.Nodes.Position={[0,0,0],[1,0,0]}';
% G3.Nodes.Position={[0,0,0],[1,0,0]}';
G12=cartesian_product(G1,G2);
G=cartesian_product(G12,G3);
figure(1)
plot(G)
L1=laplacian(G1);
L=laplacian(G);
t_final=4;

x1_0=[0;1];
y1_0=[0;0];
% z1_0=[0;0];
% x2_0=[0;0];
% y2_0=[0;1];
% z2_0=[0;0];
% x3_0=[0;0];
% y3_0=[0;0];
% z3_0=[0;1];
figure(2)
draw(L1,x1_0,y1_0,'k')
% Berechnung des dynamischen Ablaufs von einem einfachen Graph mit 2 Knoten
[Tode_x1,Uode_x1]=ode45(@(t,u)ode_func(t,u,L1),[0,t_final],x1_0);
% [Tode_y1,Uode_y1]=ode45(@(t,u)ode_func(t,u,L1),[0,t_final],y1_0);
figure(3)
plot(Tode_x1,Uode_x1)
title('Berechnung des dynamischen Ablaufs von einem einfachen Graph mit 2 Knoten')
%Berechnung des Ablaufs von einem Kompositgraph
%aus Kronecker Product vom Ergebnis des einfachen Graphs
x=zeros(8,length(Uode_x1));
% y=zeros(8,length(Uode_y1));
for i=1:length(Uode_x1)
    x(:,i)=kron(kron(Uode_x1(i,:),Uode_x1(i,:)),Uode_x1(i,:));
end
% for i=1:length(Uode_y1)
%     y(:,i)=kron(kron(Uode_y1(i,:),Uode_y1(i,:)),Uode_y1(i,:));
% end
figure(4)
plot(Tode_x1,x)
title('Berechnung des Ablaufs von einem Kompositgraph aus Kronecker Product vom Ergebnis des einfachen Graphs')
% Berechnung des dynamischen Ablaufs von einem Kompositgraph aus
% Cartesian_product dreier einfachen Graphen
x_0=kron(kron(x1_0,x1_0),x1_0);
[Tode_x,Uode_x]=ode45(@(t,u)ode_func(t,u,L),[0,t_final],x_0);
% [Tode_y1,Uode_y1]=ode45(@(t,u)ode_func(t,u,L1),[0,t_final],y1_0);
figure(5)
plot(Tode_x,Uode_x)
title('Berechnung des dynamischen Ablaufs von einem Kompositgraph aus Cartesian_product dreier einfachen Graphen')
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