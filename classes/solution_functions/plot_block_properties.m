function red_sols = plot_block_properties( block_sols)
%  [xi_opt_func, eps_ips, frob_displacement_gradient, frob_green_lagrange, gl_delta_hps, gl_rot]
%     1 - 2        3              4                           5              6-7           8  <-- indices!

red_sols = [];
delta_hps_a = [];
delta_hps_b = [];
delta_hps_c = [];

for i = 1: length(block_sols.array)

    a = block_sols.array(i).shape_vec_opt;
    b = block_sols.array(i).gl_opt;
    c = block_sols.array(i).disg_opt;   
    
%     if abs( b(3) ) < 0.1
%         block_sols.array(i).lath_solution_pair
%     end
   
    if   a(3) < 0.22   &&   a(4) < 0.21   &&  a(5) < 0.42
        xx_a(i) = a(1) / a(2);
        eps_ips_a(i) = a(3);
        frob_disg_a(i) = a(4);
        frob_gl_a(i) = a(5);
        delta_hps_a = cat(1,delta_hps_a, a(6:7) );
        rot_ang_a(i) = a(8);
    end
    
    if  b(3) < 0.22   &&   b(4) < 0.21   &&  b(5) < 0.42      
        xx_b(i) = b(1) / b(2);
        eps_ips_b(i) = b(3);
        frob_disg_b(i) = b(4);
        frob_gl_b(i) = b(5);
        delta_hps_b = cat(1,delta_hps_b, b(6:7) );
        rot_ang_b(i) = b(8);
    end
       
    if c(3) < 0.22  &&  c(4) < 0.21  &&  c(5) < 0.42
        xx_c(i) = c(1) / c(2);
        eps_ips_c(i) = c(3);
        frob_disg_c(i) = c(4);
        frob_gl_c(i) = c(5);
        delta_hps_c = cat(1,delta_hps_c, c(6:7) );
        rot_ang_c (i) = c(8);
        
        red_sols = cat( 1, red_sols, block_sols.array(i) );
    end
    
end
%%
figure;
plot( xx_b, frob_gl_b,'or' );
title('ratio vs frobeniusnorm of green lagrange')

figure;
plot( xx_c, frob_disg_c,'og' );
title('ratio vs frobeniusnorm of displacement gradient')
%axis([0.4 1.8 0 inf])

figure;
plot( xx_a, eps_ips_a,'ob' );
title('ratio vs shape strain')

figure;
plot( rot_ang_b, frob_gl_b,'or');
title('rotation angle vs frobenius norm of green lagrange')

figure;
plot( rot_ang_c, frob_disg_c,'og');
title('rotation angle vs displacement gradient ')

figure;
plot( rot_ang_a, eps_ips_a,'ob');
title('rotation angle vs shape strain')

%%
figure;
subplot(2,2,1) 
plot( xx_a, frob_gl_a,'ob' );
subplot(2,2,2) 
plot( xx_b, frob_gl_b,'or' );
subplot(2,2,3) 
plot( xx_c, frob_gl_c,'og' );
title('ratio vs frobeniusnorm of green lagrange')

figure;
subplot(2,2,1) 
plot( xx_a, frob_disg_a,'ob' );
subplot(2,2,2) 
plot( xx_b, frob_disg_b,'or' );
subplot(2,2,3) 
plot( xx_c, frob_disg_c,'og' );
title('ratio vs frobeniusnorm of displacement gradient')

figure;
subplot(2,2,1) 
plot( xx_a, eps_ips_a,'ob' );
subplot(2,2,2) 
plot( xx_b, eps_ips_b,'or' );
subplot(2,2,3) 
plot( xx_c, eps_ips_c,'og' );
title('ratio vs shape strain')

figure;
subplot(2,2,1) 
plot( rot_ang_a, frob_gl_a,'ob');
subplot(2,2,2) 
plot( rot_ang_b, frob_gl_b,'or');
subplot(2,2,3) 
plot( rot_ang_c, frob_gl_c,'og');
title('rotation angle vs frobenius norm of green lagrange')

figure;
subplot(2,2,1) 
plot( rot_ang_a, eps_ips_a,'ob');
subplot(2,2,2)
plot( rot_ang_b, eps_ips_b,'or');
subplot(2,2,3) 
plot( rot_ang_c, eps_ips_c,'og');
title('rotation angle vs shape strain')

figure;
subplot(2,2,1) 
plot( rot_ang_a, frob_disg_a,'ob');
subplot(2,2,2)
plot( rot_ang_b, frob_disg_b,'or');
subplot(2,2,3) 
plot( rot_ang_c, frob_disg_c,'og');
title('rotatoin angle vs displacement gradient ')

figure;
subplot(2,2,1) 
histogram( delta_hps_a );
subplot(2,2,2) 
histogram( delta_hps_b );
subplot(2,2,3) 
histogram( delta_hps_c );

end