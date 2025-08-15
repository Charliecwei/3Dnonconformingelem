% %%%%%% 验证3D四面体上的
%%%%%P3+{la1^2la3la4 la2^2la4la1 la3^2la1la2 la4^2la2la3 po} 多项式  
%%%%%其中po=la1*la2^3;
%%%%4个面上2次多项式的积分+体上0阶矩形成的非协调元空间在四个单元(T,T1,T2,T3)下是否满足BB条件
%%%%T分别与T1,T2,T3相交于T1的1，2，3面

clc;clear;
format rat
%format shortG


%% 顶点: 1,2,3,4,5,6,7

%%%%单元顺序 [1,2,3,4]
%%%%        [7,2,4,3]
%%%%        [4,6,3,1]
%%%%        [2,1,5,4]




%%%%以T面顺序   [2,3,4]
%%%%           [3,4,1]
%%%%           [4,1,2]
%%%%           [1,2,3]



idxlamT = [1,2,3,4]; %lam on T

idxF1basisT = [1,2,3,4,5,6]; %Face basis 顺序 T
idxF1basisT1 = [1,3,2,6,5,4]; %Face basis 顺序 T1
idxlamT1 = [1,2,4,3]; %lam on T1


idxF2basisT = [1,2,3,4,5,6]; %Face basis 顺序 T
idxF2basisT2 = [1,3,2,6,5,4]; %Face basis 顺序 T2
idxlamT2 = [4,2,3,1]; %lam on T2


idxF3basisT = [1,2,3,4,5,6]; %Face basis 顺序 T
idxF3basisT3 = [1,3,2,6,5,4]; %Face basis 顺序 T3
idxlamT3 = [2,1,3,4]; %lam on T3







%% 生成F面上基函数
alpha=Duoxiangs3(3);

beta = [2,0,0;
        1,1,0;
        1,0,1;
        0,2,0;
        0,1,1;
        0,0,2];
    
delta=[0,0,0,0];
alphao=[2,0,1,1;1,2,0,1;1,1,2,0;0,1,1,2;1,3,0,0];
alpha=[alpha;alphao];
E1=shengchengMian(alpha,beta);
E2=shengchengTi(alpha,delta);
E=[E1;E2];
G=inv(E);   


%% 计算F面上基函数的体积分
delta = Duoxiangs3(1);
ET = shengchengTi(alpha,delta);
ET = (ET*G)';

for i=1:25
    for j=1:4
        if abs(ET(i,j))<1e-13
            ET(i,j)=0;
        end
    end
end
ET = round(1000000*ET)/1000000;

ETF1 = ET(1:6,:);     %F1面
ETF2 = ET(7:12,:);    %F2面
ETF3 = ET(13:18,:);   %F3面
ETF4 = ET(19:24,:);   %F4面

%%%%%%%%个面上体积分值
ETF1T = ETF1(idxF1basisT,idxlamT);
ETF1T1 = ETF1(idxF1basisT1,idxlamT1);


ETF2T = ETF2(idxF2basisT,idxlamT);
ETF2T2 = ETF2(idxF2basisT2,idxlamT2);

ETF3T = ETF3(idxF3basisT,idxlamT);
ETF3T3 = ETF3(idxF3basisT3,idxlamT3);






%%%% compute ker
ker_XF1T = null(ETF1T','r');
ker_XF1T1 = null(ETF1T1','r');

ker_XF2T = null(ETF2T','r');
ker_XF2T2 = null(ETF2T2','r');

ker_XF3T = null(ETF3T','r');
ker_XF3T3 = null(ETF3T3','r');





%%%% compute int(phi*nalba p) on T
ker_XF1T1
ker_XF1T1'*ETF1T



ker_XF2T2
ker_XF2T2'*ETF2T


ker_XF3T3
ker_XF3T3'*ETF3T




%%%%%%%%系数满足方程:S*C=0
%%%%%%%%C = [c1 c2 c3 c4 c12  c13  c23]'
% % %     S = [0   1  -1   0  -4   4   0;
% % %          0   1   0  -1  -4   0   4;
% % %          0   0  -5   5   0  12 -12;
% % %         11   0   5   0   4 -16  -4;
% % %         -9   1   0   0   8  -4   4;
% % %          0  -1   0  -7 -12   8  12];
%syms c1 c2 c3 c4 c12 c13 c23 real;
%C = [c1,c2,c3,c4,c12,c13,c23]';
%S*C     
% % 进而推出C = c*[1,1,1,1,1,1,1]';
