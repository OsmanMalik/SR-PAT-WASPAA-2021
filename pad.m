function retMat = pad(type, M, no_blocks, padding)
%pad I use this function to appropriately pad matrices when computing the
%implicit D2D forward operator and its adjoint.

block_size = size(M,1)/no_blocks;
if mod(block_size,1) ~= 0
    error('ERROR: Number of rows of M must be divisible by no_blocks')
end

retMat = zeros(size(M,1)+no_blocks*padding, size(M,2));
if type == 1
    for q = 1:no_blocks
        retMat(1 + (block_size + padding)*(q-1) : block_size*q + padding*(q-1), :) = M(1 + block_size*(q-1) : block_size*q, :);
        retMat(1 + block_size*q + padding*(q-1) : (block_size + padding)*q) = 0;
    end
elseif type == 2
    p_ceil = ceil(padding/2);
    p_floor = floor(padding/2);
    for q = 1:no_blocks
        retMat(1 + (block_size + padding)*(q-1) : (block_size + padding)*(q-1) + p_ceil, :) = 0;
        retMat(1 + (block_size + padding)*(q-1) + p_ceil : block_size*q + padding*(q-1) + p_ceil, :) = M(1 + block_size*(q-1) : block_size*q, :);
        retMat(1 + block_size*q + padding*(q-1) + p_ceil : (block_size + padding)*q, :) = 0;
    end
else
    error('ERROR: Invalid type')
end

end
