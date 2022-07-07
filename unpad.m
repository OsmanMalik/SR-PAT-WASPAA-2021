function H = unpad(type, H, no_blocks, padding)
%unpad This function is used to undo the paddning introduced by the
%function "pad". I use this function when computing the implicit D2D 
%forward operator and its adjoint.

Q = no_blocks;
P = padding;
L = size(H,1)/Q - P;
if type == 1
    for q = 1:Q
        H(1 + L*q : L*q + P, :) = [];
    end
elseif type == 2
    for q = 1:Q
        H(1+L*(q-1) : L*(q-1) + ceil(P/2), :) = [];
        H(1+L*q : L*q + floor(P/2), :) = [];
    end
end

end
