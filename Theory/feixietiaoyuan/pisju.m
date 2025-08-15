function s=pisju(i,j)
%第j个全局坐标在第i个面的位置
if j<i
    s=4+j-i;
else
   s=j-i;
end
