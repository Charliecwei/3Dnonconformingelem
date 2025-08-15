clc; clear;
%% 20220620


alphas=Duoxiangs3(3);
betaMian=Duoxiangs2(1);
index = find(sum(alphas == 1,2).*sum(alphas == 2,2));
alphas = alphas(index,:);
EF = shengchengMian(alphas,betaMian);

syms F12 F13 F14 F23 F24 F21 F34 F31 F32 F41 F42 F43 real;

F = [F12 F13 F14 F23 F24 F21 F34 F31 F32 F41 F42 F43];
x = null(EF(:,2:12)','r');

F*x

