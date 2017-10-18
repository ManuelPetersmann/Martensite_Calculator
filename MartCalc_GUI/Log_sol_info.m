function Log_sol_info(hObject, handles, red_sols,crit )
if (size( red_sols.array, 2)==1) && isempty(red_sols.array(1).F1)
    updateLog_MartCalc(hObject, handles,'No Solution fullfilling specified criteria');
    handles.lath_solutions = false;
else
    updateLog_MartCalc(hObject, handles,['Solutions reduced to : ', num2str(length(red_sols.array)), crit ] );
end
end
