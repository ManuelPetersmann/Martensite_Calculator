% If I additionally include the condition that lambda_2 of the mixing must be exactly 0 this is too restrictive!!
% also it is not necessary. Ok in the range of all solutions, ca. tolerance 1.e-3
%         reformulation of FIRST MINORS RULE - IPS constraint on the boundary:
%         Bisection to find solution for xi between 0 and 1.
        x_left  = 0.02; % linke intervallschranke -> 1-x rechte intervallschranke
        x_right = 1.- x_left;
        y2_left  = mix_y2( x_left, F1, F2);
        y2_right = mix_y2( x_right, F1, F2);
        % yy = [y2_1, y2_2];
        %if ( any( yy > 1)  &&  any( yy < 1) )
        if (y2_right -1. ) * (y2_left -1. ) < 0
            while true
                x = (x_right + x_left)/2.;
                mid = mix_y2( x, F1, F2) ;
                if abs(mid -1.) < lambda2_tol
                    %    fracs(isol) = x;
                    blocks = blocks +1;
                    break
                end
                %
                if (y2_right -1. ) * (mid -1. ) < 0
                    x_left = x;
                else
                    x_right = x;
                    y2_right = mix_y2( x_right, F1, F2);
                end
            end
        else % no local minimum for block ! Do not mix laths to blocks
            continue
        end
