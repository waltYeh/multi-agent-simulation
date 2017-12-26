function draw_graph(L,ux,uy,arg)
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
% hold off
end