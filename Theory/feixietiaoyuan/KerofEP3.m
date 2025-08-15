%% 2022.04.12
% Let E_{ijk}(u) = \frac{1}{|F_i|}\int_{F_i}u \lambda_j\lambda_k ds
% compute the E=sum_{i,j,k=1,...,4,i\neq j,k} c_{ijk}E_{ijk}
% such that E(p)=0, for all p\in P_2(K);

clear
format rat

%%
syms E122 E123 E124 E133 E134 E144 real;  
syms E233 E234 E231 E244 E241 E211 real; 
syms E344 E341 E342 E311 E312 E322 real;   
syms E411 E412 E413 E422 E423 E433 real;


alphas=Duoxiangs3(3);
betaMian=Duoxiangs2(2);
E=shengchengMian(alphas,betaMian);
G = rref(E');

idx1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20];
idx2 = [18,21,22,23,24];
C1 = G(1:19,idx2);

E = [E122 E123 E124 E133 E134 E144 E233 E234 E231 E244 E241 E211 E344 E341 E342 E311 E312 E322 E411 E412 E413 E422 E423 E433];


E1_symbol = E(idx1)';
E2_symbol = E(idx2)';


E_symbol = C1'*E1_symbol-E2_symbol

% check the rank of la1^2la3la4 la2^2la4la1 la3^2la1la2
% la4^2la2la3,la1^2la2^2;
alphao=[2,0,1,1;1,2,0,1;1,1,2,0;0,1,1,2;2,2,0,0];
E = shengchengMian(alphao,betaMian);
E1 = E(idx1,:);
E2 = E(idx2,:);
E_main = C1'*E1-E2;
for i=1:5
    for j=1:5
        if abs(E_main(i,j))<1e-14
            E_main(i,j)=0;
        end
    end
end
E_main
r = rank(E_main)