function s=pssi(i,j)
%第i个面的第j个重心坐标对应于全局位置
s=rem(i+j,4);
if s==0
    s=4;
end