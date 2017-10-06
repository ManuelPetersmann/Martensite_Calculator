function out = check_input_uitable(in)
%
read_row = 1;
out = zeros(1,3); % default return
for row = 1:size(in,1) % loop rows of table
   %in(row,:) 
   %any( isnan(in(row,:) ) )
   %any( in(row,:) )
   if any( isnan(in(row,:) ) )  ||  ~any( in(row,:) ) % check if any entry is nan (when a user delets the default zeros...) or if all elements are 0
       continue
   else
       out(read_row,:) = in(row,:);
       read_row = read_row + 1;
   end   
end

end

