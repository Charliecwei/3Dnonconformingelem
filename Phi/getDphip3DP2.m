function Dphip = getDphip3DP2(Dlambda,lambda)
    % P2 eleme 
    p = 1;
    Dphip(:,:,10) = 4*(lambda(p,3)*Dlambda(:,:,4)+lambda(p,4)*Dlambda(:,:,3));
    Dphip(:,:,1) = (4*lambda(p,1)-1).*Dlambda(:,:,1);            
    Dphip(:,:,2) = (4*lambda(p,2)-1).*Dlambda(:,:,2);            
    Dphip(:,:,3) = (4*lambda(p,3)-1).*Dlambda(:,:,3);            
    Dphip(:,:,4) = (4*lambda(p,4)-1).*Dlambda(:,:,4);
    Dphip(:,:,5) = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
    Dphip(:,:,6) = 4*(lambda(p,1)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,1));
    Dphip(:,:,7) = 4*(lambda(p,1)*Dlambda(:,:,4)+lambda(p,4)*Dlambda(:,:,1));
    Dphip(:,:,8) = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
    Dphip(:,:,9) = 4*(lambda(p,2)*Dlambda(:,:,4)+lambda(p,4)*Dlambda(:,:,2));
end