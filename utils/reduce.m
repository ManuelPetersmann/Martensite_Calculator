    function[out] = reduce(in)
    weiter = 1;
    out = zeros(3);
    
    for i = 1:size(in,3)
       u=0;
       def1 = in(:,:,i);
       for j = 1:size(out,3)                     
            def2 = out(:,:,j);
            bool= equal(def1,def2);
            if equal(def1, def2)
              break
            end
            u = u +1 ;
       end        
        if u == size(out, 3)       
            out(:,:,weiter) = def1; %in sort schreiben
           % DeformationsgradientenSort(:,:,(size(DeformationsgradientenSort,3)+1) ) = zeros(3); %eine weitere matrix hinzuf�gen
            weiter = weiter +1;
        end
    end %�u�ere for ende 
    end % function reduzieren ende