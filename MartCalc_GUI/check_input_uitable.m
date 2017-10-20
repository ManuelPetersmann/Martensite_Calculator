function out = check_input_uitable(in)
%
read_row = 1;
out = zeros(1,3); % default return
for row = 1:size(in,1) % loop rows of table
    %in(row,:)
    %any( isnan(in(row,:) ) )
    %any( in(row,:) )
    %
    %if ~isnumeric(in(row,:))
    %    continue
    %end
    % check if all entries are integers with mod([vec])==0 PET 20.10.17 (previously isnan)
    %(when a user delets the default zeros...) or if all elements are 0
    if ( ~all(mod(in(row,:),1)== 0)  ||  ~any( in(row,:) )  )
        continue
    else
        out(read_row,:) = in(row,:);
        read_row = read_row + 1;
    end
end

end

