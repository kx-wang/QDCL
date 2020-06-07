function [ a,c,d,b ] = adjustIndices( a,c,d,b, numstates )
% input the indices, then adjust them to the lowest numbers possible that
% are >0
    
    for i = 1:20
        if min([a b c d ]) > numstates
            a = a-numstates; c = c-numstates; b = b-numstates; d = d-numstates;
        elseif min([a b c d ]) < 1
            a = a+numstates; c = c+numstates; b = b+numstates; d = d+numstates;
        else
            break;
        end 
    end

end

