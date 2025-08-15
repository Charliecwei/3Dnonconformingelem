function S=Duoxiangs2(n)
%二维n次多项式的可能组和
p=(n+2)*(n+1)/2;
S=zeros(p,3);
po=1;
for i=1:n+1
    no=i-1;
    So=Duoxiangs1(no);
    pl=size(So);
    pl=pl(1);
    S(po:po+pl-1,1)=n-no+zeros(pl,1);
    S(po:po+pl-1,2:3)=So;
    po=po+pl;
end

    