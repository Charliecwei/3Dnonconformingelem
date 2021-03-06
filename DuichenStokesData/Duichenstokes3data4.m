function pde=Duichenstokes3data4
%%%3D的stokes方程 u=curl[y.^6.*(1-y).^6*x.^5.*(1-x).^5.*z^3.*(1-z).^4
%%%，y.^5.*(1-y).^5*x.^6.*(1-x).^6.*z^3.*(1-z).^4 , 0];
%%% p=0;
%%% Dirichlet边界条件


exactDu = struct('u1',@u1,'u2',@u2,'u3',@u3);
pde = struct('f', @f, 'exactp', @exactp, 'exactu',@exactu,'exactDu',exactDu,'g_D',@exactu,'g_N',@g_N);

%%子函数

    function s=f(p)
       x=p(:,1); y=p(:,2);z=p(:,3);
       
      s = [ - (y - 1/2).*(z - 1) - 6.*x.^6.*y.^5.*(x - 1).^6.*(y - 1).^5.*(z - 1).^4 - 72.*x.^6.*y.^5.*z.*(x - 1).^6.*(y - 1).^5.*(z - 1).^3 - 12.*x.^6.*y.^5.*z.^3.*(2.*z - 2).*(x - 1).^6.*(y - 1).^5 - 90.*x.^4.*y.^5.*z.^2.*(x - 1).^6.*(y - 1).^5.*(z - 1).^4 - 120.*x.^4.*y.^5.*z.^3.*(x - 1).^6.*(y - 1).^5.*(z - 1).^3 - 216.*x.^5.*y.^5.*z.^2.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 - 288.*x.^5.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^3 - 60.*x.^6.*y.^3.*z.^2.*(x - 1).^6.*(y - 1).^5.*(z - 1).^4 - 80.*x.^6.*y.^3.*z.^3.*(x - 1).^6.*(y - 1).^5.*(z - 1).^3 - 150.*x.^6.*y.^4.*z.^2.*(x - 1).^6.*(y - 1).^4.*(z - 1).^4 - 200.*x.^6.*y.^4.*z.^3.*(x - 1).^6.*(y - 1).^4.*(z - 1).^3 - 90.*x.^6.*y.^5.*z.^2.*(x - 1).^4.*(y - 1).^5.*(z - 1).^4 - 60.*x.^6.*y.^5.*z.^2.*(x - 1).^6.*(y - 1).^3.*(z - 1).^4 - 108.*x.^6.*y.^5.*z.^2.*(x - 1).^6.*(y - 1).^5.*(z - 1).^2 - 120.*x.^6.*y.^5.*z.^3.*(x - 1).^4.*(y - 1).^5.*(z - 1).^3 - 80.*x.^6.*y.^5.*z.^3.*(x - 1).^6.*(y - 1).^3.*(z - 1).^3, 6.*x.^5.*y.^6.*(x - 1).^5.*(y - 1).^6.*(z - 1).^4 - (x - 1/2).*(z - 1) + 72.*x.^5.*y.^6.*z.*(x - 1).^5.*(y - 1).^6.*(z - 1).^3 + 12.*x.^5.*y.^6.*z.^3.*(2.*z - 2).*(x - 1).^5.*(y - 1).^6 + 60.*x.^3.*y.^6.*z.^2.*(x - 1).^5.*(y - 1).^6.*(z - 1).^4 + 80.*x.^3.*y.^6.*z.^3.*(x - 1).^5.*(y - 1).^6.*(z - 1).^3 + 150.*x.^4.*y.^6.*z.^2.*(x - 1).^4.*(y - 1).^6.*(z - 1).^4 + 200.*x.^4.*y.^6.*z.^3.*(x - 1).^4.*(y - 1).^6.*(z - 1).^3 + 90.*x.^5.*y.^4.*z.^2.*(x - 1).^5.*(y - 1).^6.*(z - 1).^4 + 120.*x.^5.*y.^4.*z.^3.*(x - 1).^5.*(y - 1).^6.*(z - 1).^3 + 216.*x.^5.*y.^5.*z.^2.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 + 288.*x.^5.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^3 + 60.*x.^5.*y.^6.*z.^2.*(x - 1).^3.*(y - 1).^6.*(z - 1).^4 + 90.*x.^5.*y.^6.*z.^2.*(x - 1).^5.*(y - 1).^4.*(z - 1).^4 + 108.*x.^5.*y.^6.*z.^2.*(x - 1).^5.*(y - 1).^6.*(z - 1).^2 + 80.*x.^5.*y.^6.*z.^3.*(x - 1).^3.*(y - 1).^6.*(z - 1).^3 + 120.*x.^5.*y.^6.*z.^3.*(x - 1).^5.*(y - 1).^4.*(z - 1).^3, 36.*x.^5.*y.^5.*z.*(x - 1).^6.*(y - 1).^5.*(z - 1).^4 - 36.*x.^5.*y.^5.*z.*(x - 1).^5.*(y - 1).^6.*(z - 1).^4 - (x - 1/2).*(y - 1/2) - 36.*x.^5.*y.^6.*z.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 + 36.*x.^6.*y.^5.*z.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 - 120.*x.^3.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^6.*(z - 1).^4 + 120.*x.^3.*y.^5.*z.^3.*(x - 1).^6.*(y - 1).^5.*(z - 1).^4 - 120.*x.^3.*y.^6.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 - 300.*x.^4.*y.^5.*z.^3.*(x - 1).^4.*(y - 1).^6.*(z - 1).^4 + 540.*x.^4.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 - 300.*x.^4.*y.^6.*z.^3.*(x - 1).^4.*(y - 1).^5.*(z - 1).^4 - 120.*x.^5.*y.^3.*z.^3.*(x - 1).^5.*(y - 1).^6.*(z - 1).^4 + 120.*x.^5.*y.^3.*z.^3.*(x - 1).^6.*(y - 1).^5.*(z - 1).^4 - 540.*x.^5.*y.^4.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 + 300.*x.^5.*y.^4.*z.^3.*(x - 1).^6.*(y - 1).^4.*(z - 1).^4 - 144.*x.^5.*y.^5.*z.^2.*(x - 1).^5.*(y - 1).^6.*(z - 1).^3 + 144.*x.^5.*y.^5.*z.^2.*(x - 1).^6.*(y - 1).^5.*(z - 1).^3 - 120.*x.^5.*y.^5.*z.^3.*(x - 1).^3.*(y - 1).^6.*(z - 1).^4 + 540.*x.^5.*y.^5.*z.^3.*(x - 1).^4.*(y - 1).^5.*(z - 1).^4 - 540.*x.^5.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^4.*(z - 1).^4 - 72.*x.^5.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^6.*(z - 1).^2 + 120.*x.^5.*y.^5.*z.^3.*(x - 1).^6.*(y - 1).^3.*(z - 1).^4 + 72.*x.^5.*y.^5.*z.^3.*(x - 1).^6.*(y - 1).^5.*(z - 1).^2 - 144.*x.^5.*y.^6.*z.^2.*(x - 1).^5.*(y - 1).^5.*(z - 1).^3 - 120.*x.^5.*y.^6.*z.^3.*(x - 1).^3.*(y - 1).^5.*(z - 1).^4 - 120.*x.^5.*y.^6.*z.^3.*(x - 1).^5.*(y - 1).^3.*(z - 1).^4 - 72.*x.^5.*y.^6.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^2 + 120.*x.^6.*y.^3.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 + 300.*x.^6.*y.^4.*z.^3.*(x - 1).^5.*(y - 1).^4.*(z - 1).^4 + 144.*x.^6.*y.^5.*z.^2.*(x - 1).^5.*(y - 1).^5.*(z - 1).^3 + 120.*x.^6.*y.^5.*z.^3.*(x - 1).^3.*(y - 1).^5.*(z - 1).^4 + 120.*x.^6.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^3.*(z - 1).^4 + 72.*x.^6.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^2];
 
 
       
    end


    function s=exactp(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s=-(x - 1/2).*(y - 1/2).*(z - 1);
    end


    function s=exactu(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         s = [ 3.*x.^6.*y.^5.*z.^2.*(x - 1).^6.*(y - 1).^5.*(z - 1).^4 + 4.*x.^6.*y.^5.*z.^3.*(x - 1).^6.*(y - 1).^5.*(z - 1).^3, - 3.*x.^5.*y.^6.*z.^2.*(x - 1).^5.*(y - 1).^6.*(z - 1).^4 - 4.*x.^5.*y.^6.*z.^3.*(x - 1).^5.*(y - 1).^6.*(z - 1).^3, 6.*x.^5.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^6.*(z - 1).^4 - 6.*x.^5.*y.^5.*z.^3.*(x - 1).^6.*(y - 1).^5.*(z - 1).^4 + 6.*x.^5.*y.^6.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 - 6.*x.^6.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4];
 
         
    end


    function s=u1(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s = [                                                                                                                                              18.*x.^5.*y.^5.*z.^2.*(x - 1).^6.*(y - 1).^5.*(z - 1).^4 + 24.*x.^5.*y.^5.*z.^3.*(x - 1).^6.*(y - 1).^5.*(z - 1).^3 + 18.*x.^6.*y.^5.*z.^2.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 + 24.*x.^6.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^3,                                                                                                                                              15.*x.^6.*y.^4.*z.^2.*(x - 1).^6.*(y - 1).^5.*(z - 1).^4 + 20.*x.^6.*y.^4.*z.^3.*(x - 1).^6.*(y - 1).^5.*(z - 1).^3 + 15.*x.^6.*y.^5.*z.^2.*(x - 1).^6.*(y - 1).^4.*(z - 1).^4 + 20.*x.^6.*y.^5.*z.^3.*(x - 1).^6.*(y - 1).^4.*(z - 1).^3,                                                                                                                                                                                                                                               6.*x.^6.*y.^5.*z.*(x - 1).^6.*(y - 1).^5.*(z - 1).^4 + 24.*x.^6.*y.^5.*z.^2.*(x - 1).^6.*(y - 1).^5.*(z - 1).^3 + 12.*x.^6.*y.^5.*z.^3.*(x - 1).^6.*(y - 1).^5.*(z - 1).^2];
    end



    function s=u2(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
       s = [                                                                                                                                            - 15.*x.^4.*y.^6.*z.^2.*(x - 1).^5.*(y - 1).^6.*(z - 1).^4 - 20.*x.^4.*y.^6.*z.^3.*(x - 1).^5.*(y - 1).^6.*(z - 1).^3 - 15.*x.^5.*y.^6.*z.^2.*(x - 1).^4.*(y - 1).^6.*(z - 1).^4 - 20.*x.^5.*y.^6.*z.^3.*(x - 1).^4.*(y - 1).^6.*(z - 1).^3,                                                                                                                                            - 18.*x.^5.*y.^5.*z.^2.*(x - 1).^5.*(y - 1).^6.*(z - 1).^4 - 24.*x.^5.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^6.*(z - 1).^3 - 18.*x.^5.*y.^6.*z.^2.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 - 24.*x.^5.*y.^6.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^3,                                                                                                                                                                                                                                             - 6.*x.^5.*y.^6.*z.*(x - 1).^5.*(y - 1).^6.*(z - 1).^4 - 24.*x.^5.*y.^6.*z.^2.*(x - 1).^5.*(y - 1).^6.*(z - 1).^3 - 12.*x.^5.*y.^6.*z.^3.*(x - 1).^5.*(y - 1).^6.*(z - 1).^2];
 
    end

    function s=u3(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s = [ 30.*x.^4.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^6.*(z - 1).^4 - 30.*x.^4.*y.^5.*z.^3.*(x - 1).^6.*(y - 1).^5.*(z - 1).^4 + 30.*x.^4.*y.^6.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 + 30.*x.^5.*y.^5.*z.^3.*(x - 1).^4.*(y - 1).^6.*(z - 1).^4 - 72.*x.^5.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 + 30.*x.^5.*y.^6.*z.^3.*(x - 1).^4.*(y - 1).^5.*(z - 1).^4 - 30.*x.^6.*y.^5.*z.^3.*(x - 1).^4.*(y - 1).^5.*(z - 1).^4, 30.*x.^5.*y.^4.*z.^3.*(x - 1).^5.*(y - 1).^6.*(z - 1).^4 - 30.*x.^5.*y.^4.*z.^3.*(x - 1).^6.*(y - 1).^5.*(z - 1).^4 + 72.*x.^5.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 - 30.*x.^5.*y.^5.*z.^3.*(x - 1).^6.*(y - 1).^4.*(z - 1).^4 + 30.*x.^5.*y.^6.*z.^3.*(x - 1).^5.*(y - 1).^4.*(z - 1).^4 - 30.*x.^6.*y.^4.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 - 30.*x.^6.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^4.*(z - 1).^4, 18.*x.^5.*y.^5.*z.^2.*(x - 1).^5.*(y - 1).^6.*(z - 1).^4 - 18.*x.^5.*y.^5.*z.^2.*(x - 1).^6.*(y - 1).^5.*(z - 1).^4 + 24.*x.^5.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^6.*(z - 1).^3 - 24.*x.^5.*y.^5.*z.^3.*(x - 1).^6.*(y - 1).^5.*(z - 1).^3 + 18.*x.^5.*y.^6.*z.^2.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 + 24.*x.^5.*y.^6.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^3 - 18.*x.^6.*y.^5.*z.^2.*(x - 1).^5.*(y - 1).^5.*(z - 1).^4 - 24.*x.^6.*y.^5.*z.^3.*(x - 1).^5.*(y - 1).^5.*(z - 1).^3];
    end



  function s = g_N(p)
        x=p(:,1); y=p(:,2);z=p(:,3);
        
       
        uD(:,:,1) = u1(p);
        uD(:,:,2) = u2(p);
        uD(:,:,3) = u3(p);
        
        uD1 = zeros(length(x),3);
        uD2 = zeros(length(x),3);
        uD3 = zeros(length(x),3);
        
        for i = 1:3
            uD1(:,i) = uD(:,i,1)+uD(:,1,i);
            uD2(:,i) = uD(:,i,2)+uD(:,2,i);
            uD3(:,i) = uD(:,i,3)+uD(:,3,i);
        end
        
        
        up = exactp(p);
        
        s = zeros( size(p,1),3);
        
        gg = (abs(x)<=1e-15);
        s(gg,1) =-uD1(gg,1)+up(gg);
        s(gg,2) =-uD2(gg,1);
        s(gg,3) =-uD3(gg,1);
        
        gg = (abs(1-x)<=1e-15);
        s(gg,1) = uD1(gg,1)-up(gg);
        s(gg,2) = uD2(gg,1);
        s(gg,3) = uD3(gg,1);
        
        
        
        gg = (abs(y)<=1e-15);
        s(gg,1) = -uD1(gg,2);
        s(gg,2) = -uD2(gg,2)+up(gg);
        s(gg,3) = -uD3(gg,2);
        
        gg = (abs(1-y)<=1e-15);
        s(gg,1) = uD1(gg,2);
        s(gg,2) = uD2(gg,2)-up(gg);
        s(gg,3) = uD3(gg,2);
        
        
        
        
        gg = (abs(z)<=1e-15);
        s(gg,1) = -uD1(gg,3);
        s(gg,2) = -uD2(gg,3);
        s(gg,3) = -uD3(gg,3)+up(gg);
        
        gg = (abs(1-z)<=1e-15);
        s(gg,1) = uD1(gg,3);
        s(gg,2) = uD2(gg,3);
        s(gg,3) = uD3(gg,3)-up(gg);
        
    end
        


end
         