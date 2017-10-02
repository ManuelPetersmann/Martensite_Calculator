%% calculate possible solutions and store solution objects in an object array


updateLog_MartCalc(hObject, handles, 'Determination of solutions started.')
switch calc_mech
    case 1
        %% integrated file: maraging_block_sym_doubleshear.m;
        updateLog_MartCalc(hObject, handles, 'maraging_block_sym_doubleshear - run')
        % highly symmetric mirror planes from bcc
        % {001} family
        sort_out_negatives = true;
        ms = all_from_family_perms( [0 0 1], sort_out_negatives );
        % {011} family
        ms = cat(1, ms, all_from_family_perms( [0 1 1], sort_out_negatives ) );
        martensite.mirror_planes = ms;
        %block_symmetric_doubleshear(B, cp, ms, ns, ds )
        %
    case 2
        updateLog_MartCalc(hObject, handles, 'maraging_block_sym_doubleshear_specific_slipsys_test - run')
        maraging_block_sym_doubleshear_specific_slipsys_test;
    case 3
        updateLog_MartCalc(hObject, handles, 'maraging_MarescaCurtin_test - run')
        maraging_MarescaCurtin_test;
    case 4
        updateLog_MartCalc(hObject, handles, 'maraging_multiple_shears - run')
        maraging_multiple_shears;
    case 5
        updateLog_MartCalc(hObject, handles, 'maraging_variable_doubleshear - run')
        maraging_variable_doubleshear;
end



switch calculation_input
    % case 0 - cannot happen - error already occurs in function 'independent_slipsystems'
    case 1
martensite.IPS_solutions = doubleshear_variable_shear_mags( martensite.U, austenite.slip_planes, ...
                         austenite.slip_directions, martensite.cp, martensite.slip_planes, martensite.slip_directions);
    case 2                     
martensite.IPS_solutions = doubleshear_variable_shear_mags( martensite.U, austenite.slip_planes, austenite.slip_directions);
    case 3
martensite.IPS_solutions = doubleshear_variable_shear_mags( martensite.U, martensite.slip_planes, ...
                                                           martensite.slip_directions, martensite.cp);
end