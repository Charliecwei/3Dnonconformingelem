clear
format rat
alpha = Duoxiangs2(3);
betaMian = Duoxiangs1(2);
betaTi = Duoxiangs2(0);
EMian = GeneraldDshengchengMian(alpha,betaMian);
ETi = GeneraldDshengchengTi(alpha,betaTi);

E = [EMian;ETi];

G = inv(E);


for i=1:10
    for j=1:10
        if abs(G(i,j))<1e-14
            G(i,j)=0;
        end
    end
end


%%
%%2019.6.12
%%计算基函数g,及其导数 Dg
 syms la1 la2 la3  real;
s = [la1 la2 la3];
f = prod(power(s,alpha),2);
f = f';
g = f*G;
g = g';



syms Da1 Da2 Da3 real;
S = [Da1 Da2 Da3];
Df = zeros(10,1);
for j = 1:3
    alphas = alpha;
    alphas(:,j) = alphas(:,j)-1;
    Df = Df+alpha(:,j).*prod(power(s,alphas),2).*S(j);
end
 Df = Df';
 Dg = Df*G;
 Dg = Dg';
