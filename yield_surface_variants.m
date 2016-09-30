
for sig_1 = 0. : 0.1 : 100
    for sig_2 = 0. : 0.1 : 100
        for sig_2 = 0. : 0.1 : 100
            sig = [ sig_1 0. 0.; 0. sig_2 0.; 0. 0. sig_3];
            rss = trace(sig' * epsilon(:,:,i) )
        end
    end
end