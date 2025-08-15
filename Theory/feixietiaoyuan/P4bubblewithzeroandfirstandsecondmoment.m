clear
clc
format rat

%%% P4
%%% 面上0,1,2阶矩自由度下的bubble function
%%% 2022.05.30





%%
betaMian=Duoxiangs2(2);
alpha=Duoxiangs3(3);

%% la1^4
alpha1 = [4,0,0,0];

%%
alpha2 = [3,1,0,0
          1,3,0,0
          1,1,2,0
          1,1,0,2];

% alpha2 = [2,2,0,0];  

% alpha2 = [2,2,0,0
%           2,0,2,0
%           2,0,0,2];
      

% alpha2 = [1,3,0,0
%           1,0,3,0
%           1,0,0,3];


% alpha2 = [2,1,1,0
%           2,1,0,1
%           2,0,1,1];



alphas = [alpha;alpha1];      
E1=shengchengMian(alphas,betaMian); 
ker_x = null(E1,'r');
g1 = Jihanhsu(alphas,ker_x);



alphas = [alpha;alpha2];      
E2=shengchengMian(alphas,betaMian); 
ker_x = null(E2,'r');
ker_x = rref(ker_x')';
g2 = Jihanhsu(alphas,ker_x);


E3=shengchengMian(alpha,betaMian); 
ker_x = null(E3,'r');
ker_x = rref(ker_x')';
g3 = Jihanhsu(alpha,ker_x);
%%
% betaMian=Duoxiangs2(3);
% alphas=Duoxiangs3(4);
% 
% E3=shengchengMian(alphas,betaMian); 
% ker_x = null(E3,'r');
% g3 = Jihanhsu(alphas,ker_x);




%%
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