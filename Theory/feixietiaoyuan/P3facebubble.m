alpha=Duoxiangs3(3);
alphao=[1,1,2,0;0,1,1,2];
alpha = [alpha;alphao];
%alpha = [3,0,0,0;2,0,0,0;1,0,0,0;0,0,0,0];
beta = Duoxiangs2(2);

E=shengchengMian(alpha,beta);

E = E(1:end-6,:);

x = null(E,'c');
x = rref(x')';

m = size(x);
for i=1:m(1)
    for j=1:m(2)
        if abs(x(i,j))<1e-10
            x(i,j)=0;
        end
    end
end

g = Jihanhsu(alpha,x);

X = [0       
       1       
       1       
      -1       
      -1       
       0       
       0       
      -1       
       0       
       1       
      -2       
       9       
      19       
       9       
    -136       
     -19       
      -2       
      19       
     -19       
       2       
       0       
     336];
 
g = Jihanhsu(alpha,X);

Dg = Grad_g(alpha,X);
syms Da1 Da2 Da3 Da4
Dg = subs(Dg,[Da1,Da2,Da3,Da4],[0,0,-1,1]);

%% check Bb condition for P1 perp
format rat
alpha = Duoxiangs3(2);
alpha1 = Duoxiangs3(1);

beta = [1,1,0,0;      %la1*la2
        1,0,1,0;      %la1*la3
        0,1,1,0;      %la2*la3
        1,1,1,0;      %la1*la2*la3
        1,0,2,0;      %la1*la3^2
        0,1,2,0;      %la2*la3^2
        0,0,0,2;      %la4^2
        0,0,0,1;      %la4
        0,0,0,0];     %1
E = shengchengTi(alpha,beta);
E1 = shengchengTi(alpha1,beta);

Es = E;
idx = logical(sum(alpha == 2,2));
Es(:,idx) = 2*Es(:,idx) - E1;
Es(:,~idx) = 4*Es(:,~idx);

P = [[12,1,0,0,12,0,0,12,1,12];
     [12,0,1,0,12,0,1,12,0,12];
     [12,0,0,1,12,1,0,12,0,12]]';

EP1 = [-1,0,1,0,0,0,0,0,0;
       -1,1,0,0,0,0,0,0,0;
       -1,0,0,0,0,0,0,0,0;
        0,0,0,-2,0,1,0,0,0;
        0,0,0,-2,1,0,0,0,0;
        0,0,0,-2,0,0,0,0,0;
        0,0,0,0,0,0,105,-90,15];

EP2 =  [-1,0,1,0,0,0,0,0,0;
        -1,1,0,0,0,0,0,0,0;
        1,0,0,0,0,0,0,0,0;
        0,0,0,-2,0,1,0,0,0;
        0,0,0,-2,1,0,0,0,0;
        0,0,0, 2,0,0,0,0,0;
        0,0,0,0,0,0,-105,90,-15];
Es1 = EP1*Es*P;
Es2 = EP2*Es*P;
EE = [Es1,Es2];





%EP3 = [-1 2 4 -2 -10 -18 136 5 -38 20 -336];
EP3 = [-2,0,2,2,10,-154,98,25,-76,25,672,-336];
EP4 = -EP3;
alpha0 = [0,1,1,1;
          0,1,0,2];
beta = [alpha;alpha0];


E = shengchengTi(alpha,beta);
E1 = shengchengTi(alpha1,beta);
Es = shengchengTi(alpha,beta);
idx = logical(sum(alpha == 2,2));
Es(:,idx) = 2*Es(:,idx) - E1;
Es(:,~idx) = 4*Es(:,~idx);
Es3 = EP3*Es*P;
Es4 = EP4*Es*P;

EE = [EE;
      Es3,Es4];
  
rank(EE) 
