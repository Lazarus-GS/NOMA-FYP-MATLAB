function out = BPSK_demod (array,n_bits)

out = zeros(1,n_bits);

    for k =1:n_bits
        if(array(k)> 0)
            out(k) = 1;
        elseif (array(k) < 0) 
            out(k) = 0;
        end
    end
    
end