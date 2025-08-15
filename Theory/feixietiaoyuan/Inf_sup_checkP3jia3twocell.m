% %%%%%% 验证3D四面体上的
%%%%%P3+{la1^2la3la4 la2^2la4la1 la3^2la1la2 la4^2la2la3 po} 多项式  
%%%%%其中po=la1^2la2^2;
%%%%4个面上2次多项式的积分+体上0阶矩形成的非协调元空间在两个单元下是否满足BB条件

clc;clear;
%format rat
format shortG


%% case 1
%%%%单元顺序 [1,2,3,4]
%%%%%%%%%%  [1,3,2,5]
%%%%面顺序   [1,2,3]
%%%%%%%% not inf sup



%% case 2
%%%%单元顺序 [1,2,3,4]
%%%%%%%%%%  [2,3,1,5]
%%%%面顺序   [1,2,3]
%%%%%%%% inf sup

%% case 3
%%%%单元顺序 [1,2,3,4]
%%%%%%%%%%  [3,1,2,5]
%%%%面顺序   [1,2,3]
%%%%%%%% inf sup

%% case 4
%%%%单元顺序 [1,2,3,4]
%%%%%%%%%%  [4,1,5,2]




Face_index1 = 4; %相交面为T1的第四的一个面
idxFbasisT1 = [1,2,3,4,5,6]; %Face basis 顺序 T1
idxlamT1 = [1,2,3,4]; %lam on T1

%%
%对应case1
%Face_index2 = 4;%相交面为T2的第四的一个面
%idxFbasisT2 = [1,3,2,6,5,4]; %Face basis 顺序 T2
%idxlamT2 = [1,3,2,4]; %lam on T2
%%

%%
%对应case2
%Face_index2 = 4;%相交面为T2的第四的一个面
%idxFbasisT2 = [4,5,2,6,3,1]; %Face basis 顺序 T2
%idxlamT2 = [2,3,1,4]; %lam on T2

%%对应case3
%Face_index2 = 4;%相交面为T2的第四的一个面
%idxFbasisT2 = [6,3,5,1,2,4];
%idxlamT2 = [3,1,2,4];


%%对应case4
Face_index2 = 3;%相交面为T2的三的一个面
idxFbasisT2 = [6,3,5,1,2,4];
idxlamT2 = [4,1,3,2];





%% 生成F4面上基函数
alpha=Duoxiangs3(3);
beta = [2,0,0;1,1,0;1,0,1;0,2,0;0,1,1;0,0,2];
delta=[0,0,0,0];
alphao=[2,0,1,1;1,2,0,1;1,1,2,0;0,1,1,2;2,2,0,0];
alpha=[alpha;alphao];
E1=shengchengMian(alpha,beta);
E2=shengchengTi(alpha,delta);
E=[E1;E2];
G=inv(E);   


%% 计算F4面上基函数的体积分
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

ET1 = ET((Face_index1-1)*6:Face_index1*6,:);
ET2 = ET((Face_index2-1)*6:Face_index2*6,:);

ET1 = ET1(idxFbasisT1,idxlamT1);
ET2 = ET2(idxFbasisT2,idxlamT2);



A = ET1';
B = ET2';
C = [A;B];


ker_XE1 = null(A,'r');
ker_XE2 = null(B,'r');
ker_XE1_XE2 = null(C,'r');

ETphi1 = ker_XE2'*ET1

ETphi2 = ker_XE1'*ET2









% 
% [ECT1,ECT1s] = array_phi_dot_gradp(ETphi1);
% [ECT2,ECT2s] = array_phi_dot_gradp(ETphi2);
% 
% x1 = null(ECT1s,'r');
% x2 = null(ECT2s,'r');
% 
% Fx1 = x1;
% Fx2 = x2;
% Fx2(idxFbasisT2,:) = x2;
% 
% F = [Fx1,-Fx2];
% 
% F(2,:) = F(2,:) - F(1,:) - F(4,:);
% F(3,:) = F(3,:) - F(1,:) - F(6,:);
% F(5,:) = F(5,:) - F(4,:) - F(6,:);
% 
% Fs = ker_XE1_XE2'*F;
% jump_x = null(Fs,'r');
% 
% F*jump_x;
% 
% 
% 
% 
% 
% 
% 
% function [ECT,ECTs] = array_phi_dot_gradp(ET)
% %% array ET by phi*gradp with 
% %%p = c1*(2*la1-1)*la1+c2*(2*la2-1)*la2+c3*(2*la3-1)*la3
% %%+c12*4*(la1*la2+la3*la4)+c13*4*(la1*la3+la2*la4)+c23*4*(la2*la3+la1*la4)+c4*(2*la4-1)*la4
% %%the array is: c1 c12 c13 c2 c23 c3 c4
% m = size(ET);
% m = m(1);
% ECT = zeros(m*4,7);
% ECTs = zeros(2*m,6);
% for i=1:m-1
%     idx = (i-1)*4;
%     idxs = (i-1)*2;
%     ECT(idx+1,[1,2,3,5]) = ET(i,[1,2,3,4]);
%     ECT(idx+2,[4,2,3,5]) = ET(i,[2,1,4,3]);
%     ECT(idx+3,[6,2,3,5]) = ET(i,[3,4,1,2]);
%     ECT(idx+4,[7,2,3,5]) = ET(i,[4,3,2,1]);
%     
%     ECTs(idxs+1,:) = ECT(idx+1,1:6)-ECT(idx+2,1:6);
%     ECTs(idxs+2,:) = ECT(idx+1,1:6)-ECT(idx+3,1:6);
% 
% end
% end
