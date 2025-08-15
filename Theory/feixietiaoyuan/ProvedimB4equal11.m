clc; clear;
%% 20230329


alphas=Duoxiangs3(4);
betaMian=Duoxiangs2(2);
index = find(sum(alphas == 3,2)+sum(alphas == 2,2).*sum(alphas == 1,2));
alphas = alphas(index,:);
EF = shengchengMian(alphas,betaMian);

syms F122 F123 F124 F133 F134 F144 real;
syms F233 F234 F231 F244 F241 F211 real;
syms F344 F341 F342 F311 F312 F322 real;
syms F411 F412 F413 F422 F423 F433 real;

F = [F122 F123 F124 F133 F134 F144 ...
     F233 F234 F231 F244 F241 F211 ...
     F344 F341 F342 F311 F312 F322 ...
     F411 F412 F413 F422 F423 F433   ];


%%the dual basis function of la1^3la2
idx = (alphas(:,1) == 3).*(alphas(:,2) == 1);
x = null(EF(:,~idx)','r');
x = 150 * x;
F*x
x'*EF(:,~~idx)
%%the duall basis function of la1^2*la2*la3
idx = (alphas(:,1) == 2).*(alphas(:,2) == 1).*(alphas(:,3) == 1);
x = null(EF(:,~idx)','r');
x = -630 * x;
F*x
x'*EF(:,~~idx)
