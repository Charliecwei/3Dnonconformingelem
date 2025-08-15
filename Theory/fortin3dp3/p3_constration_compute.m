%%类似Fortin构造3D-P2非协调元思路，尝试构造P3非协调元
syms l0 l1 l2 real; %lambda0, lambda1, lambda2
%% generate a basis of P3
phi = cell(1,10); 
%vertex
phi{1} = 9/2*l0*(l0-1/3)*(l0-2/3);
phi{2} = 9/2*l1*(l1-1/3)*(l1-2/3);
phi{3} = 9/2*l2*(l2-1/3)*(l2-2/3);
%edge 0-1
phi{4} = 27/2*l0*l1*(l0-1/3);
phi{5} = 27/2*l0*l1*(l1-1/3);
%edge 0-2
phi{6} = 27/2*l0*l2*(l0-1/3);
phi{7} = 27/2*l0*l2*(l2-1/3);
%edge 1-2
phi{8} = 27/2*l1*l2*(l1-1/3);
phi{9} = 27/2*l1*l2*(l2-1/3);
%center
phi{10} = 27*l0*l1*l2;

%% check phi is a basis of P3
alpha = genSimplexMulIndex(2,3);
alpha = alpha([1,7,10,2,4,3,6,8,9,5],:);
E = zeros(10,10);
for i = 1:10
for j = 1:10
    E(i,j) = subs(phi{j},[l0,l1,l2],alpha(i,:)/3);
end
end

%% compute constration
E = zeros(6,10);
conf = cell(6,1);
conf{1} = l0;
conf{2} = l1;
conf{3} = l2;
conf{4} = l0*l1;
conf{5} = l0*l2;
conf{6} = l1*l2;
for i = 1:6
    for j = 1:10
        f = subs(conf{i}*phi{j},l0,1-l1-l2);
        E(i,j) = 2*int(int(f,l1,[0,1-l2]),l2,[0,1]);
    end
end
E(1:3,:) = 120*E(1:3,:);
E(4:6,:) = 840*E(4:6,:);
Es = E(:,[4,5,6,7,8,9,1,2,3,10]);
rEs = rref(Es);
