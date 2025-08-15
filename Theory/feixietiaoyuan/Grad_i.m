function D = Grad_i(G,alpha,i)
m = size(alpha);
m = m(1);
n = size(G);
n = n(2);

S = alpha./sum(alpha,2);
%S = Duoxiangs3(alpha(1,1));
%S = S./sum(S,2);
D = zeros(m,n);
index = alpha(:,i) > 0;
alphas = alpha;
alphas(index,i) = alphas(index,i)-1;
for j = 1:m
    s = S(j,:);    
    D(j,:) = (alpha(:,i).*prod(power(s,alphas),2))'*G;
end

end

