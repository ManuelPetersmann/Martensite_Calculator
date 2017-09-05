function [solutions] = ave(solutions, i1, i2, prop_string, xi)
         
P1 = solutions.array(i1).(prop_string);
P2 = solutions.array(i2).(prop_string);

        P_composite = xi* P1 + (1-xi)*P2;

        eigs(P_composite)

[U1,R1] = polardecomposition(P1)

[U2,R2] = polardecomposition(P2)

F_c2 = rotation_average( R1, R2 )* ( xi* U1 + (1-xi)*U2 )
      
eigs(F_c2)
        
end