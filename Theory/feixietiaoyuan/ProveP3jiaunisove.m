
%% 证明3D四面体上的
%%%%%P3+{la1^2la3la4 la2^2la4la1 la3^2la1la2 la4^2la2la3 la1^3la2} 多项式  
%%%%% 在面上0，1，2阶矩，体上0阶矩下唯一可解
clc; clear;
format rat

alphas=Duoxiangs3(3);
betaMian=Duoxiangs2(2);
EF = shengchengMian(alphas,betaMian);

x = null(EF','r');

syms F122 F123 F124 F133 F134 F144 F233 F234 F231 F244 F241 F211 real;
syms F344 F341 F342 F311 F312 F322 F411 F412 F413 F422 F423 F433 real;

F = [F122 F123 F124 F133 F134 F144 F233 F234 F231 F244 F241 F211 F344 F341 F342 F311 F312 F322 F411 F412 F413 F422 F423 F433];

-F*x(:,1)
-2*F*x(:,2)
-F*x(:,3)
-2*F*x(:,4)
-F*x(:,5)


% X = zeros(24,5);
% for s = 1:3
%     i = rem([s+1,s+2,s+3],4);
%     idx = find(i==0);
%     i(idx) = 4;
%     j = i([3,1,2]);
%     k = [s,s,s];
%     X(:,s) = array(i,j,k);
% end
% 
% 
% X(:,4) = 2*x(:,2);
% 
% 
% i = [1,2];
% j = [3,4];
% index = Fmap(i,j,j);
% X(index,5) = 1;
% 
% j = [4,3];
% index = Fmap(i,j,j);
% X(index,5) = -1;
% 
% i = [3,4];
% j = [4,3];
% k = [1,2];
% index = Fmap(i,j,k);
% X(index,5) = 2;
% 
% k = [2,1];
% index = Fmap(i,j,k);
% X(index,5) = -2;
% 
% 
% i = [1,2];
% j = [2,1];
% k = [4,3];
% index = Fmap(i,j,k);
% X(index,5) = 4;
% 
% k = [3,4];
% index = Fmap(i,j,k);
% X(index,5) = -4;
% 




% 
% 
% function X = array(i,j,k)
% X = zeros(24,1);
% index = Fmap(i,j,j);
% X(index) = 1;
% 
% index = Fmap(j,i,i);
% X(index) = -1;
% 
% index = Fmap(j,i,k);
% X(index) = 6;
% 
% index = Fmap(i,j,k);
% X(index) = -6;
% end
% 
% 
% function l = Fmap(i,j,k)
% l = 6*(i-1);
% n = length(i);
% beta = Duoxiangs2(2);
% for s = 1:n
%     jf = pisju(i(s),j(s));
%     kf = pisju(i(s),k(s));
%   
%     beta0 = [0,0,0];
%     beta0(jf) = beta0(jf)+1;
%     beta0(kf) = beta0(kf)+1;
%     idx = (beta(:,1)==beta0(1)).*(beta(:,2)==beta0(2)).*(beta(:,3)==beta0(3));
%     l(s) = l(s)+ find(idx);
% end
% end