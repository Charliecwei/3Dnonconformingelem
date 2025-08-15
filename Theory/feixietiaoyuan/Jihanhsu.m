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