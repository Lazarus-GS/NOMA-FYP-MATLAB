function out = BPSK_mod (array,n_bits)

out = zeros(1,n_bits);

    for k =1:n_bits
        if(array(k)==1)
            out(k) = 1;
        elseif (array(k)==0) 
            out(k) = -1;
        end
    end
    
end