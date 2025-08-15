%%
% check divB_4 = P_3/{1}
% where Bk is the bubble function in P4, such that its face moment for
% order 0,1,2,.. is zero
%%
clear
clc

k = 4;
alpha=Duoxiangs3(k);
beta = Duoxiangs2(k-2);


E = shengchengMian(alpha,beta);
m = size(E);

ker_x = null(E,'c');
error = max(max(abs(E*ker_x)))/max(max(abs(E)));
disp(['matrix solve error: ', num2str(error)]);
P4 = gen_P(4);

xs = rref(ker_x')';
xs = round(120*xs)/120;
disp(['round error: ', num2str(max(max(abs(E*xs))))]);
B4 = xs'*P4




Dg1 = Grad_i(ker_x,alpha,1);
Dg2 = Grad_i(ker_x,alpha,2);
Dg3 = Grad_i(ker_x,alpha,3);
Dg4 = Grad_i(ker_x,alpha,4);

Divg = [Dg1-Dg2,Dg1-Dg3,Dg1-Dg4];


dimDivBk = rank(Divg);




