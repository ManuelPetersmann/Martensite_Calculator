function [det_eq, grad_det] = det_mix(F1, F2, x, U)
%
FC = zeros(3);
for i = 1:3
    for j = 1:3
        FC(i,j) = lin_mix(x,i,j);
    end
end

det_eq = -det(U) + FC(1,1)*FC(2,2)*FC(3,3) + FC(1,2)*FC(2,3)*FC(3,1) + FC(1,3)*FC(2,1)*FC(3,2)...
                 - FC(1,3)*FC(2,2)*FC(3,1) - FC(1,2)*FC(2,1)*FC(3,3) - FC(1,1)*FC(2,3)*FC(3,2);
             
grad_det = x^2*( ggg(1,1)*ggg(2,2)*ggg(3,3) + ggg(1,2)*ggg(2,3)*ggg(3,1) + ggg(1,3)*ggg(2,1)*ggg(3,2)...
                -ggg(1,3)*ggg(2,2)*ggg(3,1) - ggg(1,2)*ggg(2,1)*ggg(3,3) - ggg(1,1)*ggg(2,3)*ggg(3,2); 

    function mix = lin_mix(x,i,j)
        mix = x * F1(i,j) + (1.-x)*F2(i,j);
    end

    function diff = ggg(i,j)
       %  must look different!!   diff = F1(i,j) - F2(i,j)
    end

end

