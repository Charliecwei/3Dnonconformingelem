% %%%%%% 验证3D四面体上的
%%%%%P3+{la1^2la3la4 la2^2la4la1 la3^2la1la2 la4^2la2la3 po} 多项式  
%%%%%其中po=la1^2la2^2;
%%%%4个面上2次多项式的积分+体上0阶矩形成的非协调元空间在五个单元(T1,T2,T3,T4,T5)下是否满足BB条件
%%%%T1分别与T2,T3,T4相交于T1的1，2，3，4面

clc;clear;
%format rat
format shortG


%% 顶点: 1,2,3,4,5,6,7,8

%%%%单元顺序 [1,2,3,4]
%%%%        [8,2,4,3]
%%%%        [4,7,3,1]
%%%%        [2,1,6,4]
%%%%        [1,3,2,5]



%%%%以T1面顺序  [2,3,4]
%%%%           [3,4,1]
%%%%           [4,1,2]
%%%%           [1,2,3]



idxlamT1 = [1,2,3,4]; %lam on T1

idxF1basisT1 = [1,2,3,4,5,6]; %Face basis 顺序 T1
idxF1basisT2 = [1,3,2,6,5,4]; %Face basis 顺序 T2
idxlamT2 = [1,2,4,3]; %lam on T2


idxF2basisT1 = [1,2,3,4,5,6]; %Face basis 顺序 T1
idxF2basisT3 = [1,3,2,6,5,4]; %Face basis 顺序 T3
idxlamT3 = [4,2,3,1]; %lam on T3


idxF3basisT1 = [1,2,3,4,5,6]; %Face basis 顺序 T1
idxF3basisT4 = [1,3,2,6,5,4]; %Face basis 顺序 T4
idxlamT4 = [2,1,3,4]; %lam on T4


idxF4basisT1 = [1,2,3,4,5,6]; %Face basis 顺序 T1
idxF4basisT5 = [1,3,2,6,5,4]; %Face basis 顺序 T4
idxlamT5 = [1,3,2,4]; %lam on T4





%% 生成F面上基函数
alpha=Duoxiangs3(3);
beta = [2,0,0;1,1,0;1,0,1;0,2,0;0,1,1;0,0,2];
delta=[0,0,0,0];
alphao=[2,0,1,1;1,2,0,1;1,1,2,0;0,1,1,2;2,2,0,0];
alpha=[alpha;alphao];
E1=shengchengMian(alpha,beta);
E2=shengchengTi(alpha,delta);
E=[E1;E2];
G=inv(E);   


%% 计算F面上基函数的体积分
delta = Duoxiangs3(1);
ET = shengchengTi(alpha,delta);
ET = 4*(ET*G)';

for i=1:25
    for j=1:4
        if abs(ET(i,j))<1e-10
            ET(i,j)=0;
        end
    end
end
ET = round(1000*ET)/1000;

ETF1 = ET(1:6,:);     %F1面
ETF2 = ET(7:12,:);    %F2面
ETF3 = ET(13:18,:);   %F3面
ETF4 = ET(19:24,:);   %F4面

%%%%%%%%个面上体积分值
ETF1T1 = ETF1(idxF1basisT1,idxlamT1);
ETF1T2 = ETF1(idxF1basisT2,idxlamT2);


ETF2T1 = ETF2(idxF2basisT1,idxlamT1);
ETF2T3 = ETF2(idxF2basisT3,idxlamT3);

ETF3T1 = ETF3(idxF3basisT1,idxlamT1);
ETF3T4 = ETF3(idxF3basisT4,idxlamT4);


ETF4T1 = ETF4(idxF4basisT1,idxlamT1);
ETF4T5 = ETF4(idxF4basisT5,idxlamT5);



%%%% compute ker
ker_XF1T1 = null(ETF1T1','r');
ker_XF1T2 = null(ETF1T2','r');

ker_XF2T1 = null(ETF2T1','r');
ker_XF2T3 = null(ETF2T3','r');

ker_XF3T1 = null(ETF3T1','r');
ker_XF3T4 = null(ETF3T4','r');


ker_XF4T1 = null(ETF4T1','r');
ker_XF4T5 = null(ETF4T5','r');



%%%% compute int(phi*nalba p) on T1
ker_XF1T2;
ker_XF1T2'*ETF1T1;



ker_XF2T3;
ker_XF2T3'*ETF2T1;


ker_XF3T4;
ker_XF3T4'*ETF3T1;


ker_XF4T5;
ker_XF4T5'*ETF4T1;

%%%%%%%%系数满足方程:S*C=0
%%%%%%%%C = [c1 c12 c13 c2 c23  c3  c4]'
% %     S = [0  -4   4   1   0  -1   0;
% %          0  -4   0   1   4   0  -1;
% %          0   0  -4   0   4   3  -3;
% %          5   4  -8   0  -4   3   0;
% %          3   0   4  -3  -4   0   0;
% %          3   4  -4   0  -8   0   5;
% %          1   0   4  -1  -4   0   0;
% %          1   4   0   0  -4  -1   0]
% % 进而推出C = c*[1,1,1,1,1,1,1]';
