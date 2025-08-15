% %function P3jia3
% %%%%%% 验证3D四面体上的
% %%%%%P3+{la1^2la3la4 la2^2la4la1 la3^2la1la2 la4^2la2la3 po} 多项式  
% %%%%%其中po=(la1-la2)(la2-la3)(la3-la4)(la4-la1)
% %%%%4个面上2次多项式的积分下是否存在唯一
% %%%%date 2019.4.27
% 
% format rat
% 
% alpha=Duoxiangs3(3);
% %beta=Duoxiangs2(2);
% beta = [2,0,0;0,2,0;0,0,2;1,1,0;1,0,1;0,1,1];
% %beta=Duoxiangs2(2);
% delta=[0,0,0,0];
% 
% alphao=[2,0,1,1;1,2,0,1;1,1,2,0;0,1,1,2];
% 
% 
% 
% alpha=[alpha;alphao];
% 
% E1=shengchengMian(alpha,beta);
% E2=shengchengTi(alpha,delta);
% 
% E=[E1;E2];
% rank(E);
% 
% 
% alpha=Duoxiangs3(4);
% ES=shengchengMian(alpha,beta);
% 
% 
% %%F1
% S=ES(1:6,[28,26,25,29]);
% X1=S(:,1)+S(:,2)-S(:,3)-S(:,4);
% 
% 
% %%F2
% S=ES(7:12,[19,8,18,9]);
% X2=S(:,1)+S(:,2)-S(:,3)-S(:,4);
% 
% 
% 
% %%F3
% S=ES(13:18,[7,26,16,13]);
% X3=S(:,1)+S(:,2)-S(:,3)-S(:,4);
% 
% 
% 
% %%F4
% S=ES(19:24,[12,8,6,14]);
% X4=S(:,1)+S(:,2)-S(:,3)-S(:,4);
% 
% 
% 
% 
% 
% 
% X=[X1;X2;X3;X4;12/factorial(4+3)];
% 
% 
% 
% % % % % %改2019.10.28
% % % % % alpha0 = [1,1,1,1;2,1,0,1;1,2,1,0;0,1,2,1;1,0,1,2;2,0,2,0;0,2,0,2];
% % % % % alpha1 = [1,1,0,2;0,1,1,2;1,0,2,1;1,1,2,0;0,2,1,1;1,2,0,1;2,1,1,0;2,0,1,1];
% % % % % E0 = shengchengMian(alpha0,beta);
% % % % % E1 = shengchengMian(alpha1,beta);
% % % % % X0 = sum(E0,2)+E0(:,1)-sum(E1,2);
% % % % % 
% % % % % E0 = shengchengTi(alpha0,delta);
% % % % % E1 = shengchengTi(alpha1,delta);
% % % % % Xo = sum(E0,2)+E0(:,1)-sum(E1,2);
% % % % % X0 = [X0;Xo];
% % % % % 
% % % % % norm(X0-X)
% 
% 
% 
% E=[E,X];
% rank(E);
% 
% 
% 
% G=inv(E);
% 
% for i=1:25
%     for j=1:25
%         if abs(G(i,j))<1e-10
%             G(i,j)=0;
%         end
%     end
% end
% 
% 
% %%
% %%2019.6.12
% %%计算基函数g,及其导数 Dg
% alpha = [Duoxiangs3(3);alphao];
% 
%  syms la1 la2 la3 la4 po real;
% s = [la1 la2 la3 la4];
% f = prod(power(s,alpha),2);
% f = [f;po];
% f = f';
% g = f*G;
% g = g';
% 
% 
% 
% syms Da1 Da2 Da3 Da4 Dpo real;
% S = [Da1 Da2 Da3 Da4];
% Df = zeros(24,1);
% for j = 1:4
%     alphas = alpha;
%     alphas(:,j) = alphas(:,j)-1;
%     Df = Df+alpha(:,j).*prod(power(s,alphas),2).*S(j);
% end
%  Df = [Df;Dpo];
%  Df = Df';
%  Dg = Df*G;
%  Dg = Dg';
% 

%end



%%虽然可以构造，但基函数形式复杂------2019.11.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% 验证3D四面体上的
%%%%%P3+{la1^2la3la4 la2^2la4la1 la3^2la1la2 la4^2la2la3 po} 多项式  
%%%%%其中po=la1^2la2^2;
%%%%4个面上2次多项式的积分下是否存在唯一
clc;clear;
%format rat
format shortG
alpha=Duoxiangs3(3);
beta = [2,0,0;1,1,0;1,0,1;0,2,0;0,1,1;0,0,2];
delta=[0,0,0,0];
%alphao=[2,0,1,1;1,2,0,1;1,1,2,0;0,1,1,2;2,2,0,0];
alphao=[2,0,1,1;1,2,0,1;1,1,2,0;0,1,1,2;3,1,0,0];
%alphao = [2,1,1,0;0,2,1,1;1,0,2,1;1,1,0,2];


alpha=[alpha;alphao];

E1=shengchengMian(alpha,beta);
E2=shengchengTi(alpha,delta);

E=[E1;E2];
rank(E)

%alpha=[4,0,0,0;0,4,0,0;0,0,4,0;0,0,0,4];
% alpha = [2,2,0,0];
% ES=shengchengMian(alpha,beta);

% X = sum(ES,2);
% X = [X;shengchengTi(alpha,delta)];




G=inv(E);
for i=1:25
    for j=1:25
        if abs(G(i,j))<1e-10
            G(i,j)=0;
        end
    end
end

%G = round(10*G)/10;

if (max(max(E*G-eye(25))))>1e-12
    stop;
end




%%
%%2019.6.12
%%计算基函数g,及其导数 Dg
alpha = [Duoxiangs3(3);alphao];

 syms la1 la2 la3 la4 po real;
s = [la1 la2 la3 la4];
f = prod(power(s,alpha),2);
f = f';
g = f*G;
g = g';



syms Da1 Da2 Da3 Da4 Dpo real;
S = [Da1 Da2 Da3 Da4];
Df = zeros(25,1);
for j = 1:4
    alphas = alpha;
    alphas(:,j) = alphas(:,j)-1;
    Df = Df+alpha(:,j).*prod(power(s,alphas),2).*S(j);
end
 Df = Df';
 Dg = Df*G;
 Dg = Dg';

 %end
 
 
 
 
 
 
  %% check the BB condition        2020.9.16      ------ 2个宏单元不满足BB条件，3个宏单元ok
 node1 = [0,0,0;1,0,0;0,1,0;0,0,1];
 elem1 = [1,2,3,4];
 Dla1 = gradbasis3(node1,elem1);
 
 node2 = [0,0,0;-1,0,0;0,1,0;0,0,1];
 elem2 = [1,2,3,4];
 Dla2 = gradbasis3(node2,elem2);
 
 node3 = [0,0,0;1,0,0;0,-1,0;0,0,1];
 elem3 = [1,2,3,4];
 Dla3 = gradbasis3(node3,elem3);
 
 %node4 = [0,0,0;1,0,0;0,1,0;0,0,-1];
 
 
 P2 = Duoxiangs3(2);
 
 
 BBc = zeros(18,30);
 BBc2 = zeros(18,30);
 BBs = zeros(9,30);

for k = 1:3
    B1 = zeros(25,10);
    B2 = zeros(25,10);
    B3 = zeros(25,10);
     for ss = 1:10
         for j = 1:4
             alphas = alpha;
             alphas(:,j) = alphas(:,j)-1; 
             idx = alphas(:,j)>=0;
             B1(idx,ss) = B1(idx,ss)+Dla1(1,k,j)*alpha(idx,j).*prod(factorial(alphas(idx,:)+P2(ss,:)),2)./factorial(sum(alphas(idx,:)+P2(ss,:),2)+3);
             B2(idx,ss) = B2(idx,ss)+Dla2(1,k,j)*alpha(idx,j).*prod(factorial(alphas(idx,:)+P2(ss,:)),2)./factorial(sum(alphas(idx,:)+P2(ss,:),2)+3);
             B3(idx,ss) = B3(idx,ss)+Dla3(1,k,j)*alpha(idx,j).*prod(factorial(alphas(idx,:)+P2(ss,:)),2)./factorial(sum(alphas(idx,:)+P2(ss,:),2)+3);
         end
     end
     Bc1 = G'*B1;
     Bc2 = G'*B2;
     Bc3 = G'*B3;
     ks = (k-1)*6+1;
     BBc(ks:ks+5,1:10) = Bc1(7:12,:);
     BBc(ks:ks+5,11:20) = Bc2(7:12,:);
     
     BBc2(ks:ks+5,1:10) = Bc1(13:18,:);
     BBc2(ks:ks+5,21:30) = Bc3(13:18,:);
     
    ks = (k-1)*3+1;
    BBs(ks,1:10) = Bc1(25,:);
    BBs(ks+1,11:20) = Bc2(25,:);
    BBs(ks+2,21:30) = Bc2(25,:);
     
 
 
end
Bd1 = [BBc;BBc2;BBs];
BBs1 = BBs;
BBc1 = BBc;
rank(Bd1)

T0 = prod(factorial(P2),2)./(factorial(sum(P2,2)+3));

rank([Bd1;[T0',T0',T0']])

%% another way
for k = 1:3
    B1 = zeros(25,10);
    B2 = zeros(25,10);
     for ss = 1:10
         for j = 1:4
             alphas = P2(ss,:);
             alphas(:,j) = alphas(:,j)-1; 
             if alphas(:,j)>=0
                 B1(:,ss) = B1(:,ss)+Dla1(1,k,j)*P2(ss,j).*prod(factorial(alphas+alpha),2)./factorial(sum(alphas+alpha,2)+3);
                 B2(:,ss) = B2(:,ss)+Dla2(1,k,j)*P2(ss,j).*prod(factorial(alphas+alpha),2)./factorial(sum(alphas+alpha,2)+3);
             end
         end
     end
     
     Bc1 = G'*B1;
     Bc2 = G'*B2;
     ks = (k-1)*6+1;
     BBc(ks:ks+5,1:10) = Bc1(7:12,:);
     BBc(ks:ks+5,11:20) = Bc2(7:12,:);
     

     
     
    ks = (k-1)*2+1;
    BBs(ks,1:10) = Bc1(25,:);
    BBs(ks+1,11:20) = Bc2(25,:);
    
    
end

     ks = 1;

     BBc(ks,8) =  BBc(ks,8)+1/2;
     BBc(ks,18) = BBc(ks,18)-1/2;
     
     BBc(ks+1,9) =  BBc(ks+1,9)+1/2;
     BBc(ks+1,19) =  BBc(ks+1,19)-1/2;
     
     
     
     BBc(ks+2,3) =  BBc(ks+2,3)+1/2;
     BBc(ks+2,13) =  BBc(ks+2,13)-1/2;
     
     BBc(ks+3,10) =  BBc(ks+3,10)+1/2;
     BBc(ks+3,20) =  BBc(ks+3,20)-1/2;
     
     BBc(ks+4,4) =  BBc(ks+4,4)+1/2;
     BBc(ks+4,14) =  BBc(ks+4,14)-1/2;
     
      BBc(ks+5,1) =  BBc(ks+5,1)+1/2;
     BBc(ks+5,11) =  BBc(ks+5,11)-1/2;


Bd2 = [BBc;BBs];
rank(Bd2)
BBs2 = BBs;
BBc2 = BBc;





% Phio = G(:,25);
% k = 1;
% ss = 1;
% B1 = zeros(25,1);
% B2 = zeros(25,1);

% for j = 1:4
%              alphas1 = alpha;
%              alphas1(:,j) = alphas1(:,j)-1;
%              idx = alphas1(:,j)>=0;
%              alphas = P2(ss,:);
%              alphas(:,j) = alphas(:,j)-1; 
%              if alphas(:,j)>=0
%                  B2(:,ss) = B2(:,ss)+Dla1(1,k,j)* P2(ss,j).*prod(factorial(alphas+alpha),2)./factorial(sum(alphas+alpha,2)+3);
%              end
%                B1(idx,ss) = B1(idx,ss)+Dla1(1,k,j)*alpha(idx,j).*prod(factorial(alphas1(idx,:)+P2(ss,:)),2)./factorial(sum(alphas1(idx,:)+P2(ss,:),2)+3);
% end
% 
% Phio'*(B1+B2)
 
 
 %%
 