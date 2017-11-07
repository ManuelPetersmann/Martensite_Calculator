function K = compression_modulus( x, y, xy)

switch xy
    case 'nuG'
    K = 2*y*(1+x) / (3*(1-2*x));
    case 'nuE'
    case 'EG'
end


end

