function E=shengchengBian(alpha,beta)
%Ãæ»ý·Ö
n2=size(alpha,1);
no=size(beta,1);
n1=6*no;
E=zeros(n1,n2);
g = 0;
for i=1:3
    for j = i+1:4
        %%2020.11.24, ±ß 12£¬13£¬14£¬23£¬24£¬34Ë³Ðò
        for jj=1:no
               g=g+1;
            for r=1:n2
                   a=alpha(r,:);
                   as = a([i,j]);
                   if sum(as)<sum(a)
                        E(g,r)=0;
                   else
                       as=as+beta(jj,:);
                       E(g,r)=prod(factorial(as))/(factorial(sum(as)+1));
                   end
            end
        end
    end
end