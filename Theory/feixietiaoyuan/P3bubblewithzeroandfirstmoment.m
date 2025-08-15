clear
clc
format rat



%%% P3
%%% 面上一阶矩自由度下的bubble function
%%% 2022.05.30

alphas=Duoxiangs3(3);
betaMian=Duoxiangs2(1);
E1=shengchengMian(alphas,betaMian);

betaT=Duoxiangs3(1);
E2=shengchengTi(alphas,betaT);

betaDian = [1,0,0,0
            0,1,0,0
            0,0,1,0
            0,0,0,1];
E3 = zeros(4,20);
for i=1:4
    E3(i,:) = prod(power(betaDian(i,:),alphas),2);
end

E = [E1;E2;E3];
G=inv(E);

g = Jihanhsu(alphas,G);




%%%%%
betaMian=Duoxiangs2(1);
betaTi = Duoxiangs3(1);
alpha=Duoxiangs3(2);
alpha1 = [3,0,0,0
          0,3,0,0
          0,0,3,0
          0,0,0,3];
alpha2 = [1,1,1,0
          1,1,0,1
          1,0,1,1
          0,1,1,1];
      
alphas = [alpha;alpha1];      
E1=shengchengMian(alphas,betaMian); 
ET1 = shengchengTi(alphas,betaTi);

ker_x = 5*null(E1,'r');
T1 = (ET1*ker_x)';
T1 = T1(2:5,:);
g1 = Jihanhsu(alphas,ker_x);

alphas = [alpha;alpha2]; 
E2 = shengchengMian(alphas,betaMian);
ET2 = shengchengTi(alphas,betaTi);

ker_x = 30*null(E2,'r');
T2 = (ET2*ker_x)';
T2 = T2(2:5,:);
g2 = Jihanhsu(alphas,ker_x);

alphas = [alpha;alpha1;alpha2];
E3=shengchengMian(alphas,betaMian); 


alphas = Duoxiangs3(3);
E4=shengchengMian(alphas,betaMian); 



N =    [1 0 0 0 1 1 1 0
        0 1 0 0 1 1 0 1
        0 0 1 0 1 0 1 1
        0 0 0 1 0 1 1 1];
T = [-1/140 1/420 1/420 1/420 -1/84 -1/84 -1/84 -1/70
        1/420 -1/140 1/420 1/420 -1/84 -1/84 -1/70 -1/84
        1/420 1/420 -1/140 1/420 -1/84 -1/70 -1/84 -1/84
        1/420 1/420 1/420 -1/140 -1/70 -1/84 -1/84 -1/84];
Es = [N;T];



function g = Jihanhsu(alphas,G)
    m = size(G);
    for i=1:m(1)
        for j=1:m(2)
            if abs(G(i,j))<1e-10
                G(i,j)=0;
            end
        end
    end

    syms la1 la2 la3 la4 real;
    s = [la1 la2 la3 la4];
    f = prod(power(s,alphas),2);
    f = f';
    g = f*G;
    g=g';
end