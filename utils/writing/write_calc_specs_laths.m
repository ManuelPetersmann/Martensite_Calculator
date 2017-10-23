function write_calc_specs_laths(filename, permission, mart, reduced_solutions)
% write the calculation options: method, selection criteria, variable
% sorted for

%FileID = fid
fid = fopen(filename,permission);

fprintf(fid,'\n %s \n\n',...
'####################  CALCULATION OPTIONS FOR LATHS  ####################');

% calculation method
fprintf(fid,'%s \n \t %s \n','Calculation method: ', mart.IPS_solutions.calculation_method ); 

% User specified Selection criteria
fprintf(fid,'\n %s \n','User specified Selection criteria:');
%
if nargin < 4 % no selection criteria applied --> no handles.reduced_solutions
    fprintf(fid,'\t %s \n','No selection criteria specified.');
else % reduced solutions available
    keys = reduced_solutions.selection_criteria.keys;
    values = reduced_solutions.selection_criteria.values;
    for i = 1:length( keys )
        if strcmp(keys{i},'stepwidth')
            uneq = '>';
        else
            uneq = '<';
        end
        fprintf(fid,['\t %s ',uneq,' %3.4f \n'], keys{i}, values{i} );
    end
end


fclose(fid);

end


