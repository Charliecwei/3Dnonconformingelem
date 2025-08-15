% Try to construct the high order non-comforming elem from paper 3D CR
% family

%% comput the orighal polynomial for n order from from paper 3D CR p10
format rat;
syms x_1 x_2 x_3 real;
syms l1 l2 l3 l4 real;
% l1 = 1-(x_1+x_2+x_3), l2 = x_2, l3 = x_2, l4 = x_3;
p = 3; % order of polynimial

% Production nonconforming element, BT_nc on face 1
[BK_nc,BT_nc]=Genarl3DCR(p);

%B_nc = BT_nc;
%testBK_nc(B_nc,p)





    





function [BK_nc,BT_nc]=Genarl3DCR(p)
format rat;
syms x_1 x_2 x_3 real;
syms l1 l2 l3 l4 real;
n = p;
b = cell(n+1,1);
for k = 0:n
    b{k+1} = simplify(b_nk(n,k,x_1,x_2));
end

%check the orthonal  for n=0,1,2,3 this is right!!!
% for k = 0:n
%    Int2tringle(b{k+1},x_1,x_2)
%    Int2tringle(b{k+1}*x_1,x_1,x_2)
%    Int2tringle(b{k+1}*x_2,x_1,x_2)
%    Int2tringle(b{k+1}*x_1^2,x_1,x_2)
%    Int2tringle(b{k+1}*x_2^2,x_1,x_2)
%    Int2tringle(b{k+1}*x_1*x_2,x_1,x_2)
%end

%%  cnstruct r_p,2k
r = cell(2*floor(p/2)+1,1);

for k = 0:floor(p/2)
    rs = b{2*k+1};
    for j = 0:floor(p/2)
        M = Mp_ij(p,2*j,2*k);
        rs = rs+2*M*b{2*j+1};
    end
    r{2*k+1} = simplify(rs);
end


% % %check the orthonal  for n=0,1,2,3 this is right!!!
% % for k = 0:floor(p/2)
% %    Int2tringle(r{k+1},x_1,x_2)
% %    Int2tringle(r{k+1}*x_1,x_1,x_2)
% %    Int2tringle(r{k+1}*x_2,x_1,x_2)
% %    Int2tringle(r{k+1}*x_1^2,x_1,x_2)
% %    Int2tringle(r{k+1}*x_2^2,x_1,x_2)
% %    Int2tringle(r{k+1}*x_1*x_2,x_1,x_2)
% % end



d_trivp = floor(p/2)-floor((p-1)/3);
d_refl = floor((p+2)/3);


%% cnstruct bsym_pk
bsym_pk = cell(d_trivp,1);

if  mod(p,2)==0
    for k = 0:d_trivp-1
        bsym_pk{k+1} = r{p-2*k+1};
    end
else
     for k = 0:d_trivp-1
        bsym_pk{k+1} = r{p-2*k};
     end
end




%%  construct Lagrange baais and dual basis on Standard tetrahedron
m_alpha =  Duoxiangs3(p);
hat_Np =m_alpha/p;
hat_Np = hat_Np(:,2:end);
X = [[0,0,0,1]',[1,0,0,1]',[0,1,0,1]',[0,0,1,1]']; %Standard tetrahedron
lambda = simplify(X\([x_1,x_2,x_3,1]'))';

P = 0:p;
P = 1./factorial(P);

A = p*lambda - (0:p-1)';
A = [[1,1,1,1];A];

for i = 1:p
    A(i+1,:) = A(i,:).*A(i+1,:);
end


B = diag(P)*A;



dimp = size(hat_Np,1);
Phi = cell(dimp,1);
d = 3;
for k = 1:dimp
    phi = 1;
    for i = 1:d+1
        phi = phi*B(m_alpha(k,i)+1,i);
    end
    Phi{k} = simplify(phi);
end

% Lagrange basis
BG_pN = Phi;

%BG_pNla = subs(BG_pN,[x_1,x_2,x_3],[l2,l3,l4]);
    
% hat_T = [0,1,0;0,0,1];
% T = [1,2,1;1,1,2];
% hat_x = [[1/3;1/3],[1/2;0],[0;1/2]];


    
    % check the basis
%     XE = zeros(dimp,dimp);
%     for i = 1:dimp
%         for j = 1:dimp
%             XE(i,j) = subs(Phi{i},[x_1,x_2,x_3],hat_Np(j,:));
%         end
%     end
    




%% nonconforming bubble function on tetrahedron

vecter = [[0,0,0]',[1,0,0]',[0,1,0]',[0,0,1]'];

T = cell(4,1);
T{1} = vecter(:,[2,3,4]);
T{2} = vecter(:,[3,4,1]);
T{3} = vecter(:,[4,1,2]);
T{4} = vecter(:,[1,2,3]);

hat_T = [0,1,0;0,0,1];


BK_pknc = cell(d_trivp,1);
for ii = 0:d_trivp-1
    B_pkncs = 0;
    bsym_pknc = bsym_pk{ii+1};
    for i = 1:dimp
        k = 1;
        while k<5
            x = X_T(hat_Np(i,:)',T{k},hat_T);
            if length(x)==2
             %  subs(bsym_pknc,[x_1;x_2],x)
                B_pkncs = B_pkncs + subs(bsym_pknc,[x_1;x_2],x)*BG_pN{i};
                k = 5;
            else
                k = k+1;
            end
        end
    end
    BK_pknc{ii+1} = B_pkncs;
end



%% nonconforming bubble function on face   need to check
brefl_pk = cell(d_refl,1);
for k = 0:d_refl-1  
    bp_2k = b{2*k+1};
    brefl_pk{k+1}  =  (2*bp_2k-subs(bp_2k,[x_1,x_2],[x_2,1-x_1-x_2])-subs(bp_2k,[x_1,x_2],[1-x_1-x_2,x_1]))/3;
end


% On face 1  put vector(:,1) to origin
vecter = [[0,0,0]',[1,0,0]',[0,1,0]',[0,0,1]'];

T1 = cell(4,1);
T1{1} = vecter(:,[2,3,4]);
T1{2} = vecter(:,[1,3,4]);
T1{3} = vecter(:,[1,2,4]);
T1{4} = vecter(:,[1,2,3]);

BT_pknc = cell(d_refl,1);

for ii = 0:d_refl-1
    B_pkncs = 0;
    brefl_pknc = brefl_pk{ii+1};
    for i = 1:dimp
        k = 2;
        while k<5
            x = X_T(hat_Np(i,:)',T1{k},hat_T);
            if length(x)==2
             %  subs(bsym_pknc,[x_1;x_2],x)
                B_pkncs = B_pkncs + subs(brefl_pknc,[x_1;x_2],x)*BG_pN{i};
                k = 5;
            else
                k = k+1;
            end
        end
    end
    BT_pknc{ii+1} = B_pkncs;
end



%%

%% Construct for Barycentric coordinates 
% l1=1-x_1-x_2-x_3,l2=x_1,l3=x_2,l4=x_3;
BK_nc = BK_pknc;
BT_nc = BT_pknc;


for k = 0:d_trivp-1
    BK_nc{k+1}=subs(BK_nc{k+1},[x_1,x_2,x_3],[l2,l3,l4]);
    BK_nc{k+1}=expand(subs(BK_nc{k+1},l2+l3+l4,1-l1));
  %  f = subs(BK_nc{k+1},[l4,l1],[0,1-l2-l3])*l3*l2;
 %   Int2tringle(f,l2,l3)
end

for k = 0:d_refl-1
    BT_nc{k+1}=subs(BT_nc{k+1},[x_1,x_2,x_3],[l2,l3,l4]);
    BT_nc{k+1}=expand(subs(BT_nc{k+1},l2+l3+l4,1-l1));
 %   f = subs(BT_nc{k+1},l4,0)*l2;
    
 %    f = subs(BT_nc{k+1},[l4,l1],[0,1-l2-l3])*l3*l2;
    

 %   Int2tringle(f,l2,l3);
end


end

% Check the nonconforming bubble





%%      Some auxiliary functions

function x = X_T(x_0,T_0,T)
%% mapping hatT to T hat_T=[0,1,0;0,0,1];
% T = [T_1;T_2;T_3;];  hatx = [x_1,x_2];

T_0 = [ones(1,size(T_0,2));T_0];
x_0 = [ones(1,size(x_0,2));x_0];

    lambda = T_0\x_0;
    
    if norm(T_0*lambda-x_0)>1e2*eps
     %   fprintf('x_0 not in T_0\n');
        x = 0;
        return
    end
    
    for i = 1:3
        if abs(lambda(i))<eps
            lambda(i) = 0;
        end
    end

    x = T*lambda;

end



%%
function s = Mp_ij(p,i,j)
% the number from 3D CR family p15, 0<=i,j<=p
s = (-1)^p*pFq([-j,j+1,-i,i+1],[-p,p+2,1],1)*(2*i+1)/(p+1);
end

function p = b_nk(n,k,x_1,x_2)
% compute the origthal basis on Triangle from 3D CR family p10
p = (x_1+x_2)^k*P_nalphabeta(n-k,0,2*k+1,2*(x_1+x_2)-1)* P_nalphabeta(k,0,0,(x_1-x_2)/(x_1+x_2));
end




function p = P_nalphabeta(n,alpha,beta,x)
%compute Jacobi poliynomial from 3D CR family p10
alpha_n = factorial(alpha+n)/factorial(alpha);
%p  = alpha_n*F21(n,n+alpha+beta+1,alpha+1,(1-x)/2)/factorial(n);
a  = [-n,n+alpha+beta+1];
b = alpha+1;
z = (1-x)/2;
p = (alpha_n/factorial(n))*pFq(a,b,z);
end





function p = F21(n,b,c,z)
    % compute Gauss hypergeometric series from 3D CR family p10
    % old version can be remove
    p = 1;
    minun_k = 1;
    b_k = 1;
    c_k = 1;
    for k = 1:n
        minun_k = minun_k*(-n+k-1);
        b_k = b_k*(b+k-1);
        c_k = c_k*(c+k-1);
        p = p + (minun_k*b_k *z^k)/(c_k*factorial(k));
    end
end


function P = pFq(a,b,z)
% General hypergeometric function
k = 0;
P = 1;
Ak = 1;
Bk = 1;
k = k + 1;
Ak = Ak*prod(a+k-1);
Bk = Bk*prod(b+k-1);
while and(Ak~=0,Bk~=0)
    P = P+(Ak/Bk)*(z^k/factorial(k));
    k = k + 1;
    Ak = Ak*prod(a+k-1);
    Bk = Bk*prod(b+k-1);
end

end




function testBK_nc(B_nc,p)
syms l1 l2 l3 l4 real;
% Check nonconforming on face 4 , let la4=0, and l3=1-l1-l2;
        alpha = Duoxiangs2(p-1);
        P = [l1,l2,l3];
        P = P.^alpha;
        P = prod(P,2);
    for k = 1:length(B_nc)
        f = subs(B_nc{k},[l3,l4],[1-l1-l2,0]);
        for i = 1:length(P)
            Int2tringle(f*P(i),l1,l2)
        end
    end

end


function s = Int2tringle(f,x,y)
% compute the double int on standrand tringle {(x,y), 0<x,y<1, x+y<1}
s = int(int(f,x,0,1-y),0,1);
end
