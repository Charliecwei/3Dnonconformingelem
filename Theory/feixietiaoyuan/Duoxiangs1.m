function S=Duoxiangs1(n)
%一维n次多项式的可能组和
S=zeros(n+1,2);
for i=1:n+1
    S(i,1)=n+1-i;
    S(i,2)=i-1;
end