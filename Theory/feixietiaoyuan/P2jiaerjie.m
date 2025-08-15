%function P2jiaerjie
clear
format rat
%%% P2+{la1^3,la1^2la2,la1^2la3}
%%% 面上一阶矩连续
%%% 2019.4.27 . 失败


%%
%%%2019.4.28 .  P2+{la1 la2^2 la3,la1^2 la2,la1^2la3} ok
%%% P2+{[2,1,0,0;2,0,1,0]}+{gamma(i,:);i=12,13,14,16,18,19 中之一}
%%%%  gamma=Duoxiangs3(3);
% % %      3     0     0     0 ..............gamma(1,:)
% % %      1     1     1     0..............gamma(6,:)
% % %      1     1     0     1..............gamma(7,:)
% % %      1     0     1     1..............gamma(9,:)
% % %      0     3     0     0..............gamma(11,:)
% % %      0     1     1     1..............gamma(15,:)
% % %      0     0     3     0..............gamma(17,:)
% % %      0     0     0     3 ..............gamma(20,:)
%%%%%%%%%%%%%%%%%%%%%%属于P2向4个面的P1 的L2投影

%%% P2+{la2^2 la4,la3^2 la4}+{la2^2 la3 ,la2 la3^2}(其中之一) ok


% %%计算 P2+{la2^2 la4,la3^2 la4}+{la2^2 la3 } 下基函数
% alphao=[0,2,0,1;0,0,2,1;0,2,1,0];
% [G,g]=Jihanshu(alphao);

%%




%%
%%%%2019.4.30 
%P_2+{la2^2 la3 ,la2 la3^2}+{.....} 失败

%P_2+{la1 la2^2 , la2^2 la3}+{
%      2     0     1     0
%      2     0     0     1
%      1     0     2     0
%      1     0     0     2
%      0     0     2     1
%      0     0     1     2} ok

%P_2+{la1 la2^2 ,la1 la3^2}+{la2^2 la3,la2^2 la4,la2 la3^2,la2 la4^2,la3^2,la4,la3 la4^2}(其中之一) ok;

%计算 P_2+{la1 la2^2 ,la1 la3^2}+{la2^2 la3} 下基函数
% alphao=[1,2,0,0;1,0,2,0;0,2,1,0];
% 
% [E,r] = ranke(alphao);
% [G,g] = Jihanshu(alphao);

%%




%%
%%%%2019.5.27
% P_2+{la1^2la2,la2^2la3,la3^2la1};ok

%计算 P_2+{la1^2la2,la2^2la3,la3^2la1} 下基函数 并增加了求导运算
% alphao = [2,1,0,0;0,2,1,0;1,0,2,0];
% [E,r]=ranke(alphao);
%[G,g,Dg]=Jihanshu(alphao);

% %%%%% 2019.12.8
% % P_2+{la1la2la3,la2la3la4,la3la4la1,la4la1la2}; 失败 default.
% % P_2+{la1^2la2la4,la2^2la3la4,la3^2la4la1,la4^2la1la2} ok
% % 问题在于如何在其中找一组对称的基函数？？？？？？
% %考虑ker(E),维数为1，解 = 特解+零解的线性组合.
% alphao = [2,1,1,0;0,2,1,1;1,0,2,1;1,1,0,2];
% % X = zeros(1,14);
% % X(1,1:10) = (0.25)^2;
% % X(1,11:14) = (0.25)^4;
% [E,r] = ranke(alphao);
% Es = [E;X];
% rank(Es)
% %[G,g,Dg] = Jihanshu(alphao);
% 
%  syms la1 la2 la3 la4 real;
% 

% %%%2019.6.18 分别算在F1,F2,F3,F4上的基函数
% g1 = g([1:3,13]);
% g1 = subs(g1,la1,0);
% 
% 
% g2 = g([4:6,13]);
% g2 = subs(g2,la2,0);
% 
% 
% g3 = g([7:9,13]);
% g3 = subs(g3,la3,0);
% 
% 
% g4 = g([10:12,13]);
% g4 = subs(g4,la4,0);
% 


%2020.2.22 总结
% alphao = [3,0,0,0;0,3,0,0;0,0,3,0; 0,0,0,3];
% [E1,r]=ranke(alphao);
% E2 = fuzhizhicenter(alphao);
% E = [E1;E2];
% rank(E);




alpha=Duoxiangs3(2);
alphao1 = [2,1,0,0;2,0,1,0;2,0,0,1;0,2,1,0;0,2,0,1;0,0,2,1];
alphao2 = [1,2,0,0;1,0,2,0;1,0,0,2;0,1,2,0;0,1,0,2;0,0,1,2];
betaMian=Duoxiangs2(1);
betaTi=Duoxiangs3(1);

alphas1 = [alpha;alphao1];
alphas2 = [alpha;alphao2];
E1 = shengchengMian(alphas1,betaMian);
E2 = shengchengTi(alphas1,betaTi);
E = [E1;E2];
E1s = shengchengMian(alphas2,betaMian);
E2s = shengchengTi(alphas2,betaTi);
Es = [E1s;E2s];
rank(Es);

E = (E+Es)/2;

rank(E);


alphao = [2,1,0,0;0,2,1,0;1,0,2,0];
[E,r]=ranke(alphao);
r
%%%计算矩阵的秩,满秩13
function [E,r]=ranke(alphao)
    alphas=Duoxiangs3(2);
    betaMian=Duoxiangs2(1);
    betaTi=[0,0,0,0];
    alphas=[alphas;alphao];
    E1=shengchengMian(alphas,betaMian);
     E2=shengchengTi(alphas,betaTi);

     E=[E1;E2];
     r=rank(E);
end




%%%%%%% 计算P2+{alphao}下的基函数
function [G,g,Dg]=Jihanshu(alphao)
     betaMian=Duoxiangs2(1);
     betaTi=[0,0,0,0];

    alpha=Duoxiangs3(2);
    alpha=[alpha;alphao];


     E1=shengchengMian(alpha,betaMian);
     E2=shengchengTi(alpha,betaTi);

    E=[E1;E2];

    G=inv(E);

    for i=1:13
        for j=1:13
            if abs(G(i,j))<1e-10
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
     
    
    %%
    %%2019.5.27增加求导
    syms Da1 Da2 Da3 Da4 real;
    S = [Da1 Da2 Da3 Da4];
    Df = zeros(13,1);
    for j = 1:4
        alphas = alpha;
        alphas(:,j) = alphas(:,j)-1;
        Df = Df+alpha(:,j).*prod(power(s,alphas),2).*S(j);
    end
     Df = Df';
     
     Dg = Df*G;
     Dg = Dg';
    
end


function E = fuzhizhicenter(alphao)
% 在中心处赋值
alpha = Duoxiangs3(2);
alpha = [alpha;alphao];
a = 1/4;
E = a.^sum(alpha,2);
E = E';
end
