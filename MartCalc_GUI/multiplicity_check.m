%
if length(red_sols.array) < 1001
    red_sols.array = multiplicity_check_due_to_slip( handles.reduced_solutions );
    updateLog_MartCalc(hObject, handles,['Solutions reduced to : ', num2str(length(red_sols.array)), crit ] );
    guidata(hObject, handles);
end