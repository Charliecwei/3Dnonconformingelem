clear
clc
format rat

%%% 2022.06.04
%% P3 
%%% 面上0，1阶矩，体上0，1阶矩，四个顶点值为自由度
%%% 唯一可解

alphas=Duoxiangs3(3);
betaMian=Duoxiangs2(1);
EF = shengchengMian(alphas,betaMian);

betaTi = Duoxiangs3(1);
ET = shengchengTi(alphas,betaTi);


betaDian = [1,0,0,0
            0,1,0,0
            0,0,1,0
            0,0,0,1];
        
ED = zeros(4,20);
for i=1:4
    ED(i,:) = prod(power(betaDian(i,:),alphas),2);
end

E = [EF;ET;ED];
G=inv(E);
gP3 = Jihanhsu(alphas,G);

G_bubble = G(:,13:20);
if max(max(abs(EF*G_bubble))) > 1e-10
    print('error');
    stop
end
gP3_bubble = Jihanhsu(alphas,G_bubble);


%%
%%% lai^2*laj, lai*laj^2, 1<=i<j<=4, 在面上0，1，阶矩下唯一可解
alpha = Duoxiangs3(3);
idx = logical(sum(alpha==2,2));
alpha = alpha(idx,:);
alpha = alpha([1,4,2,5,3,6,7,9,8,10,11,12],:);
betaMian=Duoxiangs2(1);
E=shengchengMian(alpha,betaMian);
G = inv(E);
gP12 = Jihanhsu(alpha,G);



%%
%%% la
alpha = Duoxiangs3(4);
idx = logical(sum(alpha==3,2)+sum(alpha==2,2).*sum(alpha==1,2));
alpha = alpha(idx,:);
betaMian=Duoxiangs2(2);
E=shengchengMian(alpha,betaMian);
G = inv(E);
gP24 = Jihanhsu(alpha,G);

%% P4
%%% 面上0，1, 2阶矩, 体上la1*la2*la3*la4为矩，四个顶点值，六条边中点值为自由度，唯一可解
%%% 唯一可解

alphas=Duoxiangs3(4);
betaMian=Duoxiangs2(2);
EF=shengchengMian(alphas,betaMian);


betaTi = [0,0,0,0];
      
ET = shengchengTi(alphas,betaTi);


betaDian = [1,1,0,0
            1,0,1,0
            1,0,0,1
            0,1,1,0
            0,1,0,1
            0,0,1,1]/2.0;
betaDian = [1,0,0,0
            0,1,0,0
            0,0,1,0
            0,0,0,1
            betaDian
            1,1,1,1];
        
ED = zeros(10,35);
for i=1:10
    ED(i,:) = prod(power(betaDian(i,:),alphas),2);
end        
        

        



E = [EF;ED;ET];
G=inv(E);
gP4 = Jihanhsu(alphas,G);

G_bubble = G(:,25:35);

if max(max(abs(EF*G_bubble))) > 1e-10
    print('error');
    stop
end

gP4_bubble = Jihanhsu(alphas,G_bubble);




