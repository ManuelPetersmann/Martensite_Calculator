function Fc = linmix2(x, p1, p2)

if length(x) > 1
    Fc = x(1) * p1  +  (1.-x(1)) *p2;
else
    Fc = x * p1  +  (1.-x) *p2;
end

end

