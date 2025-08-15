function P = gen_P(n)
global lam;
lam = sym('la',[4,1],'real');
alpha = Duoxiangs3(n);
P = prod(power(lam,alpha'),1)';
end