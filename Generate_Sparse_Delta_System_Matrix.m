function H = Generate_Sparse_Delta_System_Matrix(Nx, Ny, Nz, d, Nt, dt, ...
    trans_pos, xyz0, beta, c, Cp, padding)
%Generate_Sparse_Delta_System_Matrix    Creates a system matrix which
%                                       essentially transforms a density p
%                                       to a temporal signal.
%
% NOTES:
%   The resulting matrix H is similar to H0 in [1], but it does not contain
%   signals of the form constant*delta' but rather just the representation
%   itself converted to time format, i.e., as if the density started moving
%   at the speed of sound toward the transducer. If delta(r) is the
%   expansion function, then the recording is 
%       p_const/R * delta(t).
%   But in this code we further discretize this so that it has an "energy
%   preserving" property. 
%
% INPUTS:
%   Nx          - The number of spherical expansion functions along the x-
%                 axis is Nx+1. So if Nx = 10, d = 1 and xyz0 = [0 0 0], 
%                 then the spherical expansion functions will be centered 
%                 at x values 0, 1, 2, ..., 10.
%   Ny          - Same as Nx, but in the y direction.
%   Nz          - Same as Nx, but in the z direction.
%   d           - The distance between grid points, i.e., d = dx = dy = dz. 
%                 We need the distance to be the same in all directions 
%                 since we are using spherical expansion functions.
%   Nt          - Number of time steps to use. The number of time points 
%                 used will be Nt+1. So if Nt = 10 and dt = 1, then the 
%                 time array will be 0, 1, 2, ..., 10, which is 11 time 
%                 points.
%   dt          - The distance between two time steps.
%   trans_pos   - The positions of all transducers. This should be an N by
%                 3 matrix, where the i-th row contains the x, y and z 
%                 coordinates of the i-th transducer.
%
% OPTIONAL INPUTS:
%   xyz0    - The (x,y,z) coordinate of the "bottom corner" of the cuboid
%             in which the initial pressure distribution will exist. For 
%             example, if the cuboid stretches over the volym 
%             [1,3]x[2,5]x[-1,1] in x-y-z space, then xyz0 would be 
%             xyz0 = [1 2 -1]. Default value is [0 0 0].
%   beta    - Thermal coefficient of volume expansion. Default value is 
%             0.000214.
%   c       - The speed of sound. Default value is 1500 m/s.
%   Cp      - The specific heat capacity of the medium at 
%             constant pressure. Default value is 1.
%   padding - The amount of padding inserted at the end of each time
%             series. In the case of multiple transducers, the padding will
%             be added after the rows of H consisting of each tranducer.
%
% REFERENCES:
% [1]   K. Wang, S. A. Ermilov, R. Su, P-H Brecht, A. A. Oraevsky, M. A.
%       Anastasio. An Imaging Model Incorporating Ultrasonic Transducer
%       Properties for Three-Dimensional Optoacoustic Tomography. IEEE
%       Trans Med Imaging 30(2), 2011.

% =========================================================================
% HANDLING FUNCTION INPUTS AND SETTINGS DEFAULT VALUES
% =========================================================================

% Also put error handling of inputs here
% Also see
% https://www.mathworks.com/matlabcentral/answers/21012-default-parameter-if-the-user-input-is-empty

mni = 7; % The number of required inputs

% Verify necessary inputs
if nargin < mni
    error('Too few inputs provided')
else
    % Check input Nx
    if ~Is_Nonneg_Integer(Nx)
        error('Nx must be a nonnegative integer.')
    end
    
    % Check input Ny
    if ~Is_Nonneg_Integer(Ny)
        error('Ny must be a nonnegative integer.')
    end
    
    % Check input Nz
    if ~Is_Nonneg_Integer(Nz)
        error('Nz must be a nonnegative integer.')
    end
    
    % Check input d
    if ~(d>0)
        error('The distance between grid points d must be positive.')
    end
    
    % Check input Nt
    if ~Is_Nonneg_Integer(Nt)
        error('Nt must be a nonnegative integer.')
    end
    
    % Check input dt
    if ~(dt>0)
        error('The distance between time points dt must be positive.')
    end
    
    % Check input trans_pos
    if isempty(trans_pos)
        error('Provide at the position of at least one transducer')
    elseif size(trans_pos,2) ~= 3
        error('Transducer position incorrectly formatted')
    end
end

% Verify and set default values for optional inputs
% Check optional input xyz0
if nargin < mni + 1 || isempty(xyz0)
    xyz0 = [0 0 0];
elseif length(xyz0) ~= 3
    error('xyz0 should be a vector of length 3')
end

% Check optional input beta
if nargin < mni + 2 || isempty(beta)
    beta = 0.000214;
elseif ~(isreal(beta) && isscalar(beta))
    error('beta must be a real scalar.')
end

% Check optional input c
if nargin < mni + 3 || isempty(c)
    c = 1500;
elseif ~(isreal(c) && isscalar(c))
    error('c must be a real scalar.')
end

% Check optional input Cp
if nargin < mni + 4 || isempty(Cp)
    Cp = 1;
elseif ~(isreal(Cp) && isscalar(Cp))
    error('Cp must be a real scalar.')
end

% Check optional input padding
if nargin < mni + 5 || isempty(padding)
    padding = 0;
end

% =========================================================================
% COMPUTE SPECIAL DELTA SYSTEM MATRIX H
% =========================================================================

Q = size(trans_pos,1); % Number of transducers
N = (Nx+1)*(Ny+1)*(Nz+1); % Number of spherical expansion functions used

t = linspace(0, Nt*dt, Nt+1)'; % Vector of all time positions

x = linspace(0, Nx*d, Nx+1) + xyz0(1);
y = linspace(0, Ny*d, Ny+1) + xyz0(2);
z = linspace(0, Nz*d, Nz+1) + xyz0(3);

p_const = beta*c^2 / (2*Cp);
h_dist = d/2; % Half distance between adjacent grid points
row = [];
col = [];
val = [];
for q = 1:Q
    for zi = 1:Nz+1
        for yi = 1:Ny+1
            for xi = 1:Nx+1
                R = sqrt((trans_pos(q,1) - x(xi)).^2 + (trans_pos(q,2) - y(yi)).^2 + (trans_pos(q,3) - z(zi)).^2);
                signal_qn = zeros(size(t));
                tstar = R/c; % The exact time at which signal should be recorded
                tidxs = find(t>=tstar); % All time indexes such that t(tidxs) >= tstar
                if ~isempty(tidxs) % Only assign values if the pulse arrives within the time span we're interested in
                    tidx = tidxs(1); % The smallest time index such that t(tidx) >= tstar
                    alpha = (t(tidx) - tstar)/dt; % Find weight to put in signal_qn(tidx-1)
                    signal_qn(tidx) = (1-alpha) * p_const/R; % Assign value
                    signal_qn(tidx-1) = alpha * p_const/R; % Assign value
                end
                nnz_idx = find(signal_qn ~= 0);
                row = [row; nnz_idx + (Nt+1+padding)*(q-1)];
                col = [col; (zi + (Nz+1)*(yi-1) + (Nz+1)*(Ny+1)*(xi-1))*ones(length(nnz_idx),1)];
                val = [val; signal_qn(nnz_idx)];
            end
        end
    end
end
H = sparse(row, col, val, (Nt+1+padding)*Q, N);

end % End of Generate_System_Matrix
