%% calculate possible solutions and store solution objects in an object array

updateLog_MartCalc(hObject, handles, 'Determination of solutions started.')
switch handles.popup_calc_mech.Value
    case 1 
        updateLog_MartCalc(hObject, handles, 'variable doubleshear incremental optimization lath level - started')
        switch martensite.considered_plasticity
            % case 0 - cannot happen - error already occurs in function 'independent_slipsystems'
            case 1
                martensite.IPS_solutions = doubleshear_variable_shear_mags( martensite.U, martensite.slip_planes, ...
                    martensite.slip_directions, martensite.cp);                
            case 2
                martensite.IPS_solutions = doubleshear_variable_shear_mags( martensite.U, austenite.slip_planes, austenite.slip_directions);
            case 3
                martensite.IPS_solutions = doubleshear_variable_shear_mags( martensite.U, austenite.slip_planes, ...
                    austenite.slip_directions, martensite.cp, martensite.slip_planes, martensite.slip_directions);
        end
    case 2
        %% integrated file: maraging_block_sym_doubleshear.m;
        updateLog_MartCalc(hObject, handles, 'direct block approach, mirrorsym. & equal double-shears - started')
        % highly symmetric mirror planes from bcc
        % {001} family
        sort_out_negatives = true;
        ms = all_from_family_perms( [0 0 1], sort_out_negatives );
        % {011} family
        ms = cat(1, ms, all_from_family_perms( [0 1 1], sort_out_negatives ) );
        martensite.mirror_planes = ms;
        switch martensite.considered_plasticity
            case 1 
                martensite.IPS_solutions = block_symmetric_doubleshear(martensite.U, ms, martensite.slip_planes, martensite.slip_directions, martensite.cp);
            case 2
                martensite.IPS_solutions = block_symmetric_doubleshear(martensite.U, ms, austenite.slip_planes, austenite.slip_directions);                                                      
            case 3
                martensite.IPS_solutions = block_symmetric_doubleshear(martensite.U, ms, martensite.slip_planes, martensite.slip_directions, martensite.cp,...
                                                                                         austenite.slip_planes, austenite.slip_directions);
        end
        %
        %% other cases could be added here
        %     case 4
        %         updateLog_MartCalc(hObject, handles, 'multiple shears incremental minimization - started')
        %         maraging_multiple_shears;
        %     case 5
        %         updateLog_MartCalc(hObject, handles, '_MarescaCurtin_test - run')
        %         maraging_MarescaCurtin_test;
end


