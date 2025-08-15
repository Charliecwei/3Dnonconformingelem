%% 2022.04.12
% Let E_{ij}(u) = \frac{1}{|F_j|}\int_{F_i}u \lambda_j ds
% compute the E=sum_{i,i=1,...,4,i\neq j} c_{ij}E_{ij}
% such that E(p)=0, for all p\in P_2(K);

clear
format rat

%%
syms E12 E13 E14 real;  
syms E23 E24 E21 real; 
syms E34 E31 E32 real;   
syms E41 E42 E43 real;


alphas=Duoxiangs3(2);
betaMian=Duoxiangs2(1);
E=shengchengMian(alphas,betaMian);
G = rref(E');

E1 = [E12,E13,E14,E23,E24,E21,E34,E31,E41]';
E2 = [E32,E42,E43]';

C1 = G(1:9,[9,11,12]);

E_symbol = C1'*E1-E2

% check the rank of \lam1^2\lam2, lam2^2lam3, lam3^2lam1
alphao = [2,1,0,0;0,2,1,0;1,0,2,0];
E = shengchengMian(alphao,betaMian);
E1 = E([1,2,3,4,5,6,7,8,10],:);
E2 = E([9,11,12],:);
E_main = C1'*E1-E2;
for i=1:3
    for j=1:3
        if abs(E_main(i,j))<1e-14
            E_main(i,j)=0;
        end
    end
end
E_main
r = rank(E_main)