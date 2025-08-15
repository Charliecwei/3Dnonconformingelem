format rat
alpha = Duoxiangs3(2);
beta = [0,0,0];
E1=shengchengMian(alpha,beta);

beta = [0,0];
E2=shengchengBian(alpha,beta);

E = [E1;E2];

lambdam1 = [1/3,1/3,1/3];

lambda = [0,1/3,1/3,1/3;1/3,0,1/3,1/3;
                   1/3,1/3,0,1/3;1/3,1/3,1/3,0;
                 1/2,1/2,0,0;1/2,0,1/2,0;1/2,0,0,1/2;
                 0,1/2,1/2,0;0,1/2,0,1/2;0,0,1/2,1/2];
             
E = shengchengDian(alpha,lambda);


rank(E)

G=inv(E);

for i=1:10
    for j=1:10
        if abs(G(i,j))<1e-16
            G(i,j)=0;
        end
    end
end


syms la1 la2 la3 la4 real;
s = [la1 la2 la3 la4];
f = prod(power(s,alpha),2);
f = f';
g = f*G;
g = g';
