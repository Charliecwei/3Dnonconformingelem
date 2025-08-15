%% check the leamma and theore in paper.
clear; clc;

%% Define the lam.
global lam la1 la2 la3 la4;
lam = sym('la',[4,1],'real');
la1 = lam(1,1);
la2 = lam(2,1);
la3 = lam(3,1);
la4 = lam(4,1);

%% Define the index of cell on face 1, face2, face3, face4
% cfi(i,:) stand for face i's vertex index
global cfi;
cfi = [2,3,4;
       3,4,1; 
       4,1,2; 
       1,2,3];
   
%% Define the 1 order face moment and 2 order face moment
global Fm1;
global Fm2;
Fm1 = cell(4,4);    % Fm_1{i,j}(f) = (\int_{F_i}f*lam_j dF_i)/area(F_i);
Fm2 = cell(4,4,4);  % Fm_2{i,j,k}(f) = (\int_{F_i}f*lam_j*lam_k dF_i)/area(F_i);
for i = 1:4
    for j = 1:4
        Fm1{i,j} = @(f)(m_1(f,i,j));
        for k = 1:4
            Fm2{i,j,k} = @(f)(m_2(f,i,j,k));
        end
    end
end

%% Define the Bubble function space: B3 and B4
global B3;
B3 = sym("b3", [8, 1],  'real');
for i = 1:4
    B3(i,1) = 5*power(lam(i,1),3) - 5*power(lam(i,1),2) +lam(i,1);
    B3(i+4,1) = -2*power(lam(i,1),2) + 4*lam(i,1) - 2 + 3*sum(power(lam(cfi(i,:),1),2)) + 30*prod(lam(cfi(i,:),1));
end

global B4;
B4 = sym("b4", [11,1],  'real');
for i = 1:4
   B4(i,1) = 14*power(lam(i,1),4) - 21*power(lam(i,1),3) + 9*power(lam(i,1),2) - lam(i,1);
end
kn = 5;
for i = 1:4
    K = setdiff([1,2,3,4],i);
    B4(kn) = 7*power(lam(i),4) - 6*power(lam(i),3) - sum(2*power(lam(K),3)*lam(i) - 9*power(lam(K),2)*power(lam(i),2));
    kn = kn+1;
end
B4(kn,1) = prod(lam);
kn = kn + 1;
i = 1;
for j = i+1:3
    K = setdiff([1,2,3,4],[i,j]);
    k = K(1);
    l = K(2);
    lai = lam(i,1);
    laj = lam(j,1);
    lak = lam(k,1);
    lal = lam(l,1);
    la = lai+laj;
    B4(kn,1) =   28*power(la, 4) - 53*power(la, 3) + 27*power(la, 2) - 18*lai*laj + 3*lak*lal ...
               + la*(21*lai*laj - 2 - 3*power(lak+lal,2) - 21*lak*lal);
    kn = kn + 1;
end

%% Define O3 and O4
global O3;
O3 = sym("O3", [12,1],  'real');
k = 1;
for i = 1:4
    for j = i+1:4
        O3(k,1) = lam(i,1)*lam(i,1)*lam(j,1);
        O3(k+1,1) = lam(i,1)*lam(j,1)*lam(j,1);
        k = k+2;
    end
end

global O4;
O4 = sym("O4", [24,1],  'real');
k = 1;
for i = 1:4
    for j = i+1:4
        O4(k,  1) = power(lam(i,1),3)*lam(j,1);
        O4(k+1,1) = power(lam(j,1),3)*lam(i,1);
        k = k+2;
    end
end
kn = 13;
for i = 1:4
    for j = i+1:4
        for k = j+1:4
            O4(kn,  1) = lam(i,1)*lam(j,1)*lam(k,1)*lam(i,1);
            O4(kn+1,1) = lam(i,1)*lam(j,1)*lam(k,1)*lam(j,1);
            O4(kn+2,1) = lam(i,1)*lam(j,1)*lam(k,1)*lam(k,1);
            kn = kn + 3;
        end
    end
end

%% Define b2 and b3
global b2 b3;
b2 = 2 - 4*sum(power(lam,2));
b3 = 11*sum(power(lam,3)) - 9*sum(power(lam,2)) + 72*sum(prod(lam(cfi')));


%% Check lemma and theorem in paper  
%Check_Lemma_1;
%Check_Lemma_4;
%Check_Theorem_2;
Check_Lemma_5;



%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Check Lemma 1.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Check_Lemma_1
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Begin Check Lemma 1.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
global Fm1 Fm2 cfi;
global B3 B4 O3 O4;
global lam;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% chech P3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%%%%%%%%%%%%%%%% check B3's 1-order face moment vash%%%%%%
disp("check B3's 1-order face moment vash.");
for s = 1:8
    for i = 1:4
        for j = cfi(i,:)
             val = Fm1{i,j}(B3(s,1));
             if abs(val) > 1e-15
                 disp("error: B3's 1-order face moment not vash!!!!!!!!!!");
             end
        end
    end
end
%%%%%%%%%%%%%%%%%% chech dim B3 >= 8 
% E = zeros(8,8);
% for i = 1:8
%     for j = 1:8
%         E(i,j) = tet_int(B3(i)*B3(j));
%     end
% end
%disp(['dim B3 >= ', num2str(rank(E))]);
E = zeros(8,8);
Lam = [1,0,0,0;
       0,1,0,0;
       0,0,1,0;
       0,0,0,1;
       0,1/3,1/3,1/3;
       1/3,0,1/3,1/3;
       1/3,1/3,0,1/3;
       1/3,1/3,1/3,0;];
for i = 1:8
    E(i,:) = subs(B3,lam,Lam(i,:)')';
end
a = sym("a_",[4,1],"real");
b = sym("b_",[4,1],"real");
c = [a;b];
E*c
f = B3'*c;
for i = 1:8
    subs(f,lam,Lam(i,:)')
end
disp(['dim B3 >= ', num2str(rank(E))]);
%%%%%%%%%%%%%%%%% check the 1-order face monet of O3 are linear independent
disp("check 1 order face monet of O3 are linear independent");
c = sym("c_", [4, 4],  'real');
v = 0;
s = 1;
for i = 1:4
    for j = i+1:4
        v = v + c(i,j)*O3(s);
        v = v + c(j,i)*O3(s+1);
        s = s + 2;
    end
end
%%%%%% mF
i = 1; j = 2; k = 3; l = 4; % here k is l, l is m in paper.
mF = @(f)(10*Fm1{j,i}(f) -  5*Fm1{i,j}(f) +  (Fm1{i,k}(f)+ Fm1{i,l}(f))  - 2*(Fm1{j,k}(f)+Fm1{j,l}(f))...
            + 7*(Fm1{k,j}(f)+Fm1{l,j}(f)) - 11*(Fm1{k,i}(f)+Fm1{l,i}(f)) + (Fm1{k,l}(f)+Fm1{l,k}(f)));
%%%%% mF(v)        
disp("0 = mF(v) = ");
mF(v)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% chech P4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% %%%%%%%%%%%%%%%%% check B4's 2-order face moment vash%%%%%%
disp("check B4's 2-order face moment vash.");
for s = 1:11
    for i = 1:4
        for j = cfi(i,:)
            for k = cfi(i,:)
                 val = Fm2{i,j,k}(B4(s,1));
                 if abs(val) > 1e-15
                     disp("error: B4's 2-order face moment not vash!!!!!!!!!!");
                 end
            end
        end
    end
end
% %%%%%%%%% chech dim B4 >= 11
% E = zeros(11,11);
% for i = 1:11
%     for j = 1:11
%         E(i,j) = tet_int(B4(i)*B4(j));
%     end
% end
% disp(['dim B4 >= ', num2str(rank(E))]);
E = zeros(11,11);
Lam = [1,0,0,0;
       0,1,0,0;
       0,0,1,0;
       0,0,0,1;
       1/2,1/2,0,0;
       1/2,0,1/2,0;
       1/2,0,0,1/2;
       0,1/2,1/2,0;
       0,1/2,0,1/2;
       0,0,1/2,1/2;
       1/4,1/4,1/4,1/4];
for i = 1:11
    E(i,:) = subs(B4,lam,Lam(i,:)')';
end
disp(['dim B4 >= ', num2str(rank(E))]);
% a = sym("a_",[4,1],"real");
% b = sym("b_",[4,1],"real");
% c = [a;b];
% E*c
% f = B3'*c;
% for i = 1:8
%     subs(f,lam,Lam(i,:)')
% end
% disp(['dim B3 >= ', num2str(rank(E))]);
%%%%%%%%%%%%%%%%% check the 2-order face monet of O4 are linear independent
disp("check 2 order face monet of O4 are linear independent.");
c = sym("c_", [4, 4],  'real');
d = sym("d_", [4, 4, 4], 'real');
v = 0;
s = 1;
for i = 1:4
    for j = i+1:4
        v = v + c(i,j)*O4(s);
        v = v + c(j,i)*O4(s+1);
        s = s + 2;
    end
end
for i = 1:4
    for j = i+1:4
        for k = j+1:4
            v = v + d(i,j,k)*O4(s);
            v = v + d(j,i,k)*O4(s+1);
            v = v + d(k,i,j)*O4(s+2);
            s = s + 3;
        end
    end
end
%%%%%% mF
i = 1; j = 2; k = 3; l = 4; % here k is l, l is m in paper.
mF = @(f)(2*(Fm2{i,j,j}(f)-Fm2{k,j,j}(f)-Fm2{l,j,j}(f)) + 18*(Fm2{i,k,l}(f)-Fm2{k,l,i}(f)-Fm2{l,i,k}(f))...
         +3*(Fm2{k,i,j}(f)+Fm2{l,i,j}(f)+Fm2{k,l,j}(f)+Fm2{l,j,k}(f)-Fm2{i,j,k}(f)-Fm2{i,j,l}(f))...
         +5*(Fm2{k,l,l}(f)+Fm2{l,k,k}(f)+Fm2{k,i,i}(f)+Fm2{l,i,i}(f)-Fm2{i,k,k}(f)-Fm2{i,l,l}(f)));
%%%%% mF(v)        
disp("0 = mF(v) = ");
mF(v)     
%%%%%% mF
i = 1; j = 2; k = 3; l = 4; % here k is l, l is m in paper.
mF = @(f)(6*(Fm2{j,k,k}(f)+Fm2{k,j,j}(f)) + 7*(3*Fm2{l,i,j}(f)+3*Fm2{l,i,k}(f)-Fm2{l,j,j}(f)-Fm2{l,k,k}(f))...
         +9*(Fm2{i,j,j}(f)+Fm2{i,k,k}(f)-Fm2{j,k,i}(f)-Fm2{j,k,l}(f)-Fm2{k,i,j}(f)-Fm2{k,l,j}(f))...
         +15*(2*Fm2{i,l,l}(f)-Fm2{j,i,i}(f)-Fm2{j,l,l}(f)-Fm2{k,i,i}(f)-Fm2{k,l,l}(f))...
         +18*Fm2{i,j,k}(f)-14*Fm2{l,j,k}(f) - 45*(Fm2{i,j,l}(f)+Fm2{i,k,l}(f)) + 54*(Fm2{j,l,i}(f)+Fm2{k,l,i}(f)));
%%%%% mF(v)
disp("0 = mF(v) = ");
mF(v)
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Check Lemma 1.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
end

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Check Lemma 4.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Check_Lemma_4
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Begin Check Lemma 3.2.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
global Fm1 Fm2;
global b2 b3;
%%%%%%%%%%%%%%% check the 1-order face moment of b2 is zero;%%%%%%%%%%%%
disp("check the bubble function b2 in paper.")
for i = 1:4
    for j = 1:4
        val = Fm1{i,j}(b2);
        if (abs(val)>1e-15)
            disp("error: The facemonet of b2 is not vaish!!!!!!!!!!!");
        end
    end
end
%%%%%%%%%%%%%%% check the 2-order face moment of b3 is zero;%%%%%%%%%%%%
disp("check the bubble function b3 in paper.")
for i = 1:4
    for j = 1:4
        for k = 1:4
            val = Fm2{i,j,k}(b3);
            if (abs(val)>1e-15)
                disp("error: The facemonet of b3 is not vaish!!!!!!!!!!!");
            end
        end
    end
end     
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Check Lemma 3.2.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Check Theorem 3.2.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Check_Theorem_2
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Begin Check Theorem 3.2.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
global Fm1 Fm2 la1 la2 la3 la4;

%%%%%%%%%%%%Check P3- is unisolve%%%%%%%%%%%%%%%%%
disp("%%%%%%%%Check P3- is unisolve%%%%%%%%");
P2 = gen_P(2);
d = sym('d',[10,1],'real');
w = d'*P2;
% P2 supplemental function
P2_supple = [la1*la1*la2; la2*la2*la3; la3*la3*la1];
c = sym('c',[3,1],'real');
v = w + c'*P2_supple;
% mF1 mF2 mF3
mF = cell(3,1);
mF{1,1} = @(f)(Fm1{1,2}(f) - Fm1{2,1}(f) + Fm1{3,1}(f) - Fm1{1,3}(f) + Fm1{2,3}(f) - Fm1{3,2}(f));
mF{2,1} = @(f)(Fm1{1,2}(f) - Fm1{2,1}(f) + Fm1{4,1}(f) - Fm1{1,4}(f) + Fm1{2,4}(f) - Fm1{4,2}(f));
mF{3,1} = @(f)(Fm1{1,3}(f) - Fm1{3,1}(f) + Fm1{4,1}(f) - Fm1{1,4}(f) + Fm1{3,4}(f) - Fm1{4,3}(f));
disp("0 = mF1(v) = ");
mF{1,1}(v)
disp("0 = mF2(v) = ");
mF{2,1}(v)
disp("0 = mF3(v) = ");
mF{3,1}(v)
% M = zeros(13,13);
% base = sym('base',[13,1],'real');
% base(1:10) = gen_P(2);
% base(11:13) = [la1*la1*la2; la2*la2*la3; la3*la3*la1];
% mF = cell(13,1);
% index = [1,2; 1,3; 1,4; 2,1; 2,3; 2,4; 3,1; 3,2; 3,4; 4,1; 4,2; 4,3];
% for i = 1:12
%     l = index(i,1);
%     s = index(i,2);
%     mF{i,1} = @(f)(Fm1{l,s}(f));
% end
% mF{13,1} = @(f)(tet_int(f));
% for i = 1:13
%     for j = 1:13
%         M(i,j) = mF{i}(base(j));
%     end
% end
% rank(M)


% %%%%%%%%%%%%Check P4- is unisolve%%%%%%%%%%%%%%%%%
disp("%%%%%%%%Check P4- is unisolve%%%%%%%%");
P3 = gen_P(3);
d = sym('d',[20,1],'real');
w = d'*P3;
% P3 supplemental function
P3_supple = [la1*la1*la3*la4;
             la2*la2*la4*la1;
             la3*la3*la1*la2;
             la4*la4*la2*la3;
             la1*la2*la2*la2];
c = sym('c',[5,1],'real');
v = w + c'*P3_supple;
% mF1 mF2 mF3 mF4, mF5
mF = cell(5,1);
mF{1,1} = @(f)(Fm2{1,2,2}(f) - Fm2{2,1,1}(f) + Fm2{3,1,1}(f) - Fm2{1,3,3}(f) + Fm2{2,3,3}(f) - Fm2{3,2,2}(f)...
          +6*(Fm2{1,3,4}(f) - Fm2{3,4,1}(f) + Fm2{2,4,1}(f) - Fm2{1,2,4}(f) + Fm2{3,4,2}(f) - Fm2{2,3,4}(f)));
      
mF{2,1} = @(f)(4*(Fm2{2,3,1}(f) - Fm2{2,3,4}(f) + Fm2{3,4,2}(f) - Fm2{3,1,2}(f)) + Fm2{2,4,4}(f) - Fm2{2,1,1}(f)...
                +Fm2{3,1,1}(f) - Fm2{3,4,4}(f) + 2*(Fm2{1,3,4}(f) - Fm2{1,2,4}(f) + Fm2{4,1,2}(f) - Fm2{4,1,3}(f)));
            
mF{3,1} = @(f)(Fm2{1,2,2}(f) - Fm2{2,1,1}(f)  + Fm2{4,1,1}(f) - Fm2{4,2,2}(f) + 2*(Fm2{1,4,4}(f) - Fm2{2,4,4}(f))...
          +6*(Fm2{1,2,3}(f) - Fm2{2,3,1}(f)  + Fm2{1,3,4}(f) - Fm2{3,4,1}(f) + Fm2{3,4,2}(f) - Fm2{2,3,4}(f))...
          +3*(Fm2{2,3,3}(f) - Fm2{1,3,3}(f)) + 12*(Fm2{2,4,1}(f)-Fm2{1,2,4}(f)));
      
mF{4,1} = @(f)(Fm2{1,4,4}(f) - Fm2{1,3,3}(f)  + Fm2{3,1,1}(f) - Fm2{2,1,1}(f) + Fm2{2,3,3}(f) - Fm2{3,4,4}(f)...
          +2*(Fm2{4,1,2}(f) - Fm2{4,2,3}(f)  + Fm2{1,3,4}(f) - Fm2{3,4,1}(f)) + 6*(Fm2{3,4,2}(f) - Fm2{1,2,4}(f))...
          +4*(Fm2{1,2,3}(f) - Fm2{3,1,2}(f)  + Fm2{2,4,1}(f) - Fm2{2,3,4}(f)));

mF{5,1} = @(f)(Fm2{4,1,1}(f) - Fm2{4,3,3}(f)  + 2*(Fm2{3,1,1}(f) - Fm2{1,3,3}(f) + Fm2{1,4,4}(f) - Fm2{3,4,4}(f))...
          +6*(Fm2{1,2,3}(f) - Fm2{3,1,2}(f)  + Fm2{1,3,4}(f) - Fm2{3,4,1}(f)) + 3*(Fm2{2,3,3}(f) - Fm2{2,1,1}(f))...
          +12*(Fm2{2,4,1}(f) - Fm2{1,2,4}(f)  + Fm2{3,4,2}(f) - Fm2{2,3,4}(f)));
disp("0 = mF1(v) = ");
mF{1,1}(v)
disp("0 = mF2(v) = ");
mF{2,1}(v)
disp("0 = mF3(v) = ");
mF{3,1}(v)
disp("0 = mF4(v) = ");
mF{4,1}(v)
disp("0 = mF5(v) = ");
mF{5,1}(v)
% M = zeros(25,25);
% base = sym('base',[25,1],'real');
% base(1:20) = gen_P(3);
% base(21:25) = [la1*la1*la3*la4;
%                la2*la2*la4*la1;
%                la3*la3*la1*la2;
%                la4*la4*la2*la3;
%                la1*la2*la2*la2];
% mF = cell(25,1); 
% index = [1,2,2; 1,2,3; 1,2,4; 1,3,3; 1,3,4; 1,4,4;
%          2,1,1; 2,1,3; 2,1,4; 2,3,3; 2,3,4; 2,4,4;
%          3,1,1; 3,1,2; 3,1,4; 3,2,2; 3,2,4; 3,4,4;
%          4,1,1; 4,1,2; 4,1,3; 4,2,2; 4,2,3; 4,3,3];
% for i = 1:24
%     l = index(i,1);
%     s = index(i,2);
%     m = index(i,3);
%     mF{i,1} = @(f)(Fm2{l,s,m}(f));
% end
% mF{25,1} = @(f)(tet_int(f));
% for i = 1:25
%     for j = 1:25
%         M(i,j) = mF{i}(base(j));
%     end
% end
% rank(M)

disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Check Theorem 3.2.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Check Lemma 5.1.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Check_Lemma_5
%%%% check divBk+1 = Pk/c, the proof of divBk+1 subsect Pk/c is simple,
%%%% problem is the: Pk/c subsect divBk+1
%%%% here main check k = 3
global lam B4;
%B3 = gen_B3;
f = sym("f", [4,1], 'real');
for i = 1:4
    f(i) = diff(B4(i),lam(i));
end
g = sym("g",[4,1], 'real'); %here g(i) stand for g_{i,i+1} in paper
kn = 5;
for i = 1:4
    j = mod(i,4) + 1;
    g(i) = diff(B4(kn),lam(i)) - diff(B4(kn),lam(j)); %%div(b*tij);
    g(i) = simplify(g(i));
    kn = kn + 1;
end
h1 = sym("h1", [4, 1], 'real'); %here h1(i) stand for h_{i,i+1,i+2} in paper
kn = 5;
for i = 1:4
    j = mod(i,4) + 1;
    k = mod(i+1,4) + 1;
    h1(i) = (diff(B4(kn),lam(i)) - diff(B4(kn),lam(j))) - (diff(B4(kn),lam(i)) - diff(B4(kn),lam(k))); %%div(b*tik);
    h1(i) = simplify(h1(i));
    kn = kn + 1;   
end
h2 = sym("h2", [4, 1], 'real'); %here h2(i) stand for h_{i,i+1,i+3} in paper
kn = 5;
for i = 1:4
    j = mod(i,4) + 1;
    k = mod(i+2,4) + 1;
    h2(i) = (diff(B4(kn),lam(i)) - diff(B4(kn),lam(j))) - (diff(B4(kn),lam(i)) - diff(B4(kn),lam(k))); %%div(b*tik);
    h2(i) = simplify(h2(i));
    kn = kn + 1;   
end
f4 = sym("f4",[3, 1], 'real'); %here f4(4) stand for f_{j,4}
for j = 1:3
    f4(j) = diff(B4(9),lam(j)) - diff(B4(9),lam(4)); %%div(b*t1i);
end
a = sym("a",[4,1],'real');
b = sym("b",[4,1],'real');
c = sym("c",[4,1],'real');
d = sym("d",[4,1],'real');
e = sym("e",[3,1],'real');
v = a'*f + b'*g + c'*h1 + d'*h2 + e'*f4;                
disp('take lai=1, laj=0 with i = 1,2,3,4, can get');
Lam1 = [1,0,0,0
        0,1,0,0
        0,0,1,0
        0,0,0,1
        1/4,1/4,1/4,1/4];
subs(v, lam, Lam1(1,:)')
subs(v, lam, Lam1(2,:)')
subs(v, lam, Lam1(3,:)')
subs(v, lam, Lam1(4,:)')
disp('take la1=la2=la3=la4=1/4 can get');
subs(v, lam, Lam1(5,:)')
% solve the five equation can get 
% a4 = -(a1+a2+a3),
% bi = -11/12*ai, i=1,2,3,
% b4 = 11/12*(a1+a2+a3);
disp('replace a4,b1-b4 by a1-a3');
v = subs(v, a(4),-(a(1)+a(2)+a(3)));
v = subs(v, b(1),-11*a(1)/12);
v = subs(v, b(2),-11*a(2)/12);
v = subs(v, b(3),-11*a(3)/12);
v = subs(v, b(4),11*(a(1)+a(2)+a(3))/12);
disp('check the value at above lam is 0');
subs(v, lam, Lam1(1,:)')
subs(v, lam, Lam1(2,:)')
subs(v, lam, Lam1(3,:)')
subs(v, lam, Lam1(4,:)')
subs(v, lam, Lam1(5,:)')

disp('take lai=2/5, laj=1/5 with i=1,2,3, can get');
Lam2 = [2/5,1/5,1/5,1/5
        1/5,2/5,1/5,1/5
        1/5,1/5,2/5,1/5];
subs(v, lam, Lam2(1,:)')
subs(v, lam, Lam2(2,:)')
subs(v, lam, Lam2(3,:)')
% solve the three equation can get 
% ei = -4ai, i = 1,2,3.
disp('replace e1-e3 by a1-a3');
v = subs(v, e(1),-4*a(1));
v = subs(v, e(2),-4*a(2));
v = subs(v, e(3),-4*a(3));
disp('check the value at above lam is 0');
subs(v, lam, Lam2(1,:)')
subs(v, lam, Lam2(2,:)')
subs(v, lam, Lam2(3,:)')

disp('take lai=2/4, laj=1/4 lak=1/4, lal=0 can get');
Lam3 =[2/4, 0,   1/4, 1/4
       2/4, 1/4,   0, 1/4
       2/4, 1/4, 1/4,   0
       1/4, 2/4,   0, 1/4
       1/4, 2/4, 1/4,   0
         0, 2/4, 1/4, 1/4
       1/4, 1/4, 2/4,   0
         0, 1/4, 2/4, 1/4
       1/4,   0, 2/4, 1/4
         0, 1/4, 1/4, 2/4
       1/4,   0, 1/4, 2/4];

for i = 1:11
    subs(v, lam, Lam3(i,:)')
    
end
% 提取系数
f = sym("f",[11,1],'real');
coef = [a(1:3);c;d];
for i = 1:11
    ftmp = coeffs(v,coef(i));
    f(i) = ftmp(2);
end
E = zeros(11,11);
for i = 1:11
        E(i,:) = subs(f, lam, Lam3(i,:)')';
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Auxialary function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% generate Polynomial function space of degree n for n = 2 or 3
function P = gen_P(n)
global lam;
if n == 2
    alpha = [2,0,0,0
             1,1,0,0
             1,0,1,0
             1,0,0,1
             0,2,0,0
             0,1,1,0
             0,1,0,1
             0,0,2,0
             0,0,1,1
             0,0,0,2];
elseif n == 3
    alpha = [3,0,0,0
             2,1,0,0
             2,0,1,0
             2,0,0,1
             1,2,0,0
             1,1,1,0
             1,1,0,1
             1,0,2,0
             1,0,1,1
             1,0,0,2
             0,3,0,0
             0,2,1,0
             0,2,0,1
             0,1,2,0
             0,1,1,1
             0,1,0,2
             0,0,3,0
             0,0,2,1
             0,0,1,2
             0,0,0,3];
end
P = prod(power(lam,alpha'),1)';
end

%% compute the 1 order face monet: (int_{F_i}f*lam_j d F_i)/area(F_i) 
function val = m_1(f, i, j)
global lam;
la = lam(j,1);
val = face_int(f*la,i);
end

%% compute the 2 order face monet: (int_{F_i}f*lam_j*lam_k d F_i)/area(F_i)
function val = m_2(f, i, j, k)
global lam;
la = lam(j,1)*lam(k,1);
val = face_int(f*la,i);
end


%% compute the (int_{F_i} f dF_i)/area(F_i) at the i face.
function val = face_int(f, i)
global lam;
global cfi;
la4 = lam(i,1);
la1 = lam(cfi(i,1),1);
la2 = lam(cfi(i,2),1);
la3 = lam(cfi(i,3),1);
f1 = subs(subs(f,la4,0),la1,1-la2-la3);

val = int(int(f1,la2,[0,1-la3]),la3,[0,1]);
val = 2*val;
end

%% compute the (int_T f dT)/vol(T)
function val = tet_int(f)
global la1 la2 la3 la4;
f1 = subs(f,la4,1-la1-la2-la3);
val = int(int(int(f1,la1,[0,1-la2-la3]),la2,[0,1-la3]),la3,[0,1]);
val = 6*val;
end