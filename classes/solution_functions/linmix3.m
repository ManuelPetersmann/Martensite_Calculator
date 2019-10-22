function Fc = linmix3(x, p1, p2, p3)

if length(x) < 3
    Fc = x(1) * p1  +  x(2) * p2 + (1. - x(1) - x(2)) * p3;
else
    Fc = x(1) * p1  +  x(2) * p2 + x(3) * p3;
end

end

