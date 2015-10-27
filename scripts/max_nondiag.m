mm = -1;
[row,col,val] = find(M);
for i = 1:size(row,1)
    if row(i) ~= col(i)
        mm = max(mm, abs(val(i)));
    end
end
mm