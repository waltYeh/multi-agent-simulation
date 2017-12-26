function T=delete_zero_row_col(B,del_row,del_col)
T=B;
if nargin==1
    del_row=1;
    del_col=1;
end
if del_row
    T (all(T == 0, 2),:) = [];
end
if del_col
    T (:,all(T == 0, 1)) = [];
end

end