%%
% check divB_k = P_{k-1}/{1}
% where Bk is the bubble function in Pk, such that its face moment for
% order 0,1,2,...k-2 is zero
%%
clear
clc

disp(['k   &   error   &   dim divBk   &   dim Pk-1/c   &   dim Pk  &   Number of Dof  &   dim bk']);
for k=2:16
    CheckDivBk_eq_Pk_1(k)
end


function CheckDivBk_eq_Pk_1(k)

alpha=Duoxiangs3(k);

beta = Duoxiangs2(k-2);

dimPk = factorial(k+3)/(factorial(3)*factorial(k));
dimPk_1 = factorial(k+2)/(factorial(3)*factorial(k-1));



E = shengchengMian(alpha,beta);
m = size(E);

ker_x = null(E,'c');
error = max(max(abs(E*ker_x)))/max(max(abs(E)));

n = size(ker_x);

Dg1 = Grad_i(ker_x,alpha,1);
Dg2 = Grad_i(ker_x,alpha,2);
Dg3 = Grad_i(ker_x,alpha,3);
Dg4 = Grad_i(ker_x,alpha,4);

Divg = [Dg1-Dg2,Dg1-Dg3,Dg1-Dg4];


dimDivBk = rank(Divg);


disp([num2str(k) ' & ' num2str(error) ' & ' num2str(dimDivBk) ' & ' num2str(dimPk_1-1) ' & ' num2str(dimPk) ' & ' num2str(m(1)) ' & ' num2str(n(2)) ' \\ ']);

end
