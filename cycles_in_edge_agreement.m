clear all
L=[1,-1,0,0,0;
    -1,3,-1,0,-1;
    0,-1,3,-1,-1;
    0,0,-1,2,-1;
    0,-1,-1,-1,3];
Delta=diag(diag(L));
A=L-Delta;
incidence=[-1,0,0,0,0,0;
    1,-1,0,0,0,1;
    0,0,0,1,1,-1;
    0,0,1,0,-1,0;
    0,1,-1,-1,0,0];
% L2 ist immer gleich der L
% L2=incidence*incidence';
Le=incidence'*incidence;
% die nicht-null Eigenwerte von Le und von L sind gleich
% eig_Le=eig(Le)
% eig_L=eig(L)
Adir=[0,1,0,0,0;
    0,0,0,0,1;
    0,1,0,0,0;
    0,0,1,0,0;
    0,0,1,1,0];
Deltadir=diag(Adir*ones(5,1));
Ldir=Deltadir-Adir;
Gdir=digraph(Adir);
G=graph(A);
figure(1)
plot(G)
figure(2)
plot(Gdir)
incidence_T=[-1,0,0,0,0,0;
    1,-1,0,0,0,0;
    0,0,0,1,0,0;
    0,0,1,0,0,0;
    0,1,-1,-1,0,0];
incidence_C=incidence-incidence_T;
incidence_T=delete_zero_row_col(incidence_T,0,1);
incidence_C=delete_zero_row_col(incidence_C,0,1);
% incidence_C=[0,0,0,0,0,0;
%     0,0,0,0,1,0;
%     0,0,0,1,-1,0;
%     0,0,0,-1,0,0;
%     0,0,0,0,0,0];
L_T=incidence_T*incidence_T';
L_C=incidence_C*incidence_C';
% L_T+L_C==L;
Le_T=incidence_T'*incidence_T;
Le_T=delete_zero_row_col(Le_T);
Le_C=incidence_C'*incidence_C;
Le_C=delete_zero_row_col(Le_C);

DTT_DC=incidence_T'*incidence_C;
DCT_DT=incidence_C'*incidence_T;
DTT_DC=delete_zero_row_col(DTT_DC);
DCT_DT=delete_zero_row_col(DCT_DT);
% (4.14) Pg.79
% [Le_T,DTT_DC;DCT_DT,Le_C]==Le;
T=incidence_T\incidence_C;
T=delete_zero_row_col(T);
R=[eye(length(T)),T];
reduced_Le=Le_T*R*R';
% nicht 100 perzent genau
% R'*Le_T*R==Le;
t_final = 4;
t_steps = 41;
t = linspace(0,t_final,t_steps);
ux=zeros(5,t_steps);
uy=zeros(5,t_steps);
ex=zeros(6,t_steps);
ey=zeros(6,t_steps);
ex_check=zeros(6,t_steps);
ey_check=zeros(6,t_steps);
ex_T=zeros(length(Le_T),t_steps);
ey_T=zeros(length(Le_T),t_steps);
ux(:,1)=[-1.2,-0.4,0.1,0.9,0.6]';
uy(:,1)=[-1.1,-0.4,0.6,0.9,0.1]';
ex_check(:,1)=incidence'*ux(:,1);
ey_check(:,1)=incidence'*uy(:,1);
ex_T(:,1)=incidence_T'*ux(:,1);
ey_T(:,1)=incidence_T'*uy(:,1);
ex_C=T'*ex_T(:,1);
ey_C=T'*ey_T(:,1);
ex(:,1) = [ex_T(:,1);ex_C];
ey(:,1) = [ey_T(:,1);ey_C];
figure
draw_graph(L,ux(:,1),uy(:,1),'k')
for i=1:length(t)
    T_next = t(i);
    clear Uode_x
    clear Uode_y
    if i~=1
        [Tode_x,Eode_x]=ode45(@(t,u)agreement_protocol_func(t,u,reduced_Le),[0,T_next],ex_T(:,1));
        [Tode_y,Eode_y]=ode45(@(t,u)agreement_protocol_func(t,u,reduced_Le),[0,T_next],ey_T(:,1));
        ex_T(:,i) = Eode_x(end,:);
        ey_T(:,i) = Eode_y(end,:);
        ex_C=T'*ex_T(:,i);
        ey_C=T'*ey_T(:,i);
        ex(:,i) = [ex_T(:,i);ex_C];
        ey(:,i) = [ey_T(:,i);ey_C];
        
        [Tode_x,Eode_x_check]=ode45(@(t,u)agreement_protocol_func(t,u,Le),[0,T_next],ex_check(:,1));
        [Tode_y,Eode_y_check]=ode45(@(t,u)agreement_protocol_func(t,u,Le),[0,T_next],ey_check(:,1));
        ex_check(:,i) = Eode_x_check(end,:);
        ey_check(:,i) = Eode_y_check(end,:);
    end
end
figure
plot(ex','DisplayName','ex')
figure
plot(ey','DisplayName','ey')
figure
plot(ex_check','DisplayName','ex_c')
figure
plot(ey_check','DisplayName','ey_c')
