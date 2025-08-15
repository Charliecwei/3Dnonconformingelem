function P3o2jie
%%%% 验证3D四面体上的P3 多项式 在4个顶点,4个体上的积分
%%%%4个面上1次多项式的积分下是否存在唯一
%%%%date 2019.4.27

alpha=Duoxiangs3(3);
beta=Duoxiangs2(1);
delta=Duoxiangs3(1);
 lambda=[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];

format rat



E1=shengchengDian(alpha,lambda);
E2=shengchengTi(alpha,delta);
E3=shengchengMian(alpha,beta);
E=[E1;E2;E3];
rank(E)
G=inv(E);

for i=1:20
    for j=1:20
        if abs(G(i,j))<1e-10
            G(i,j)=0;
        end
    end
end





end







