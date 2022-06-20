function states = bits_to_states(z)
% convert bitstring to set of corresponding states
    states = {};
    for c = 1:length(z):
        if z(c) == 1
            states{c} = [1;0];
        else 
            states{c} = [0;1];
        end
    end
end
        
negR = z(2*p + 3 - i);
        negL = z(2*p + 2 - i);