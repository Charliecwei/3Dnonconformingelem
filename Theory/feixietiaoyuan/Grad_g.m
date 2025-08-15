function Dg = Grad_g(alpha,G)

syms la1 la2 la3 la4 po real;
s = [la1 la2 la3 la4];

syms Da1 Da2 Da3 Da4 Dpo real;
S = [Da1 Da2 Da3 Da4];

m = size(alpha);
m = m(1);
Df = zeros(m,1);

for j = 1:4
    alphas = alpha;
    alphas(:,j) = alphas(:,j)-1;
    Df = Df+alpha(:,j).*prod(power(s,alphas),2).*S(j);
end
 Df = Df';
 Dg = Df*G;
 Dg = Dg';


end