format rat

syms l1 l2 real;
p = 5;
alpha = 0:p;
P = l1.^alpha;
beta = (0:p-1)';
as = alpha+beta;

E = 2*factorial(as)./factorial(as+2);

Es = rref(E);

x = [Es(:,end);-1];

ps = P*x;

%int(int(ps*l2,l2,0,1-l1),l1,0,1)