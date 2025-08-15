function S=Duoxiangs3(n)
%三维n次多项式的可能组和
p=(n+3)*(n+2)*(n+1)/6;
S=zeros(p,4);
po=1;
for i=1:n+1
    no=i-1;
    So=Duoxiangs2(no);
    pl=size(So);
    pl=pl(1);
    S(po:po+pl-1,1)=n-no+zeros(pl,1);
    S(po:po+pl-1,2:4)=So;
    po=po+pl;
end
