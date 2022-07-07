% This script is for doing the 1st and 2nd order reconstructions in the
% paper. The kWave reconstructions are done separetely in the following two
% scripts: 
% - kWave_reconstruction_circular.m
% - kWave_reconstruction_square.m
%
% This script needs to be run prior to running the kWave reconstruction
% scripts since those scripts use some of the outputs from this script. See
% the README file for further details.
%
% In our paper, we show results for the following parameter settings in
% this script: 
% - r = 3, 6, 12
% - transducer_setting = 'circular', 'surface array'
%
% In order to recreate the plots in Figures 4 and 5 via
% plot_waspaa_experiments.m, please first run this script for all six
% combinations of those three r values and two transducer_setting
% options. See the README file for further details.

%% Settings

plotting = false;
saving = true;
    
Nx_gen = 101;  % No. points along x-axis (generation) 
Ny_gen = 101;  % No. points along y-axis (generation) 
Nx_rec = 81;  % No. points along x-axis (reconstruction) 
Ny_rec = 81;  % No. points along y-axis (reconstruction) 
no_spec = 1000;  % Number of speckles
noise_level = .01;  % Amount of Gaussian noise added to recordings
r = 3;  % Controls speckle correlation (higher r => finer speckles)
d_gen = 4e-6*40/(Nx_gen-1); % Spatial grid spacing [m]  Default 4e-6
d_rec = 4e-6*40/(Nx_rec-1);  % Spatial grid spacing [m]  Default 4e-6
Nt = 200; % Control number of time points
dt = 1e-09; % Temporal grid spacing [s]  (4-10 times smaller than 1/(f0+FWHM/2))
xyz0 = [0 0 0]; % Positioning of grid
N_gen = Nx_gen*Ny_gen;
N_rec = Nx_rec*Ny_rec;

% Transducer settings
no_trans = 64;  % Number of transducers
transducer_spacing = 15e-6;
transducer_setting = 'surface array';  % Transducer geometry
f0 = 50e6; 
FWHM = 25e6; 
eir_length = Nt; % Length of the sampled EIR
d0 = 1500/f0;

%% Compute covariance matrix average intensity

C_E = Compute_Covariance_Matrix(Nx_rec, Ny_rec, r);

%% Create density

Rho = create_star_density(Nx_gen, Ny_gen);
Rho = Rho/max(Rho(:));

%% Create transducers

if strcmp(transducer_setting,'circular')
    rad = Nx_gen-1;
    cent = (Nx_gen-1)/2;
    trans_pos = d_gen*[cent+rad*cos(2*pi*[0:no_trans-1]'/no_trans) ... % Position of transducers
        cent+rad*sin(2*pi*[0:no_trans-1]'/no_trans) ...
        zeros(no_trans,1)]; 

elseif strcmp(transducer_setting,'line array')
    trans_pos = [linspace(xyz0(1),xyz0(1)+Nx*d,no_trans)' repmat(xyz0(2)-6*d,no_trans,1) zeros(no_trans,1)];

elseif strcmp(transducer_setting,'single')
    trans_pos = [xyz0(1)+Nx/2*d xyz0(2)-1*d 0];

elseif strcmp(transducer_setting,'surface array')  
    trans_pos = zeros(no_trans,3);
    no_trans_sqrt = sqrt(no_trans);
    x_mid = d_gen*(Nx_gen-1)/2;
    y_mid = d_gen*(Ny_gen-1)/2;
    surf_side_mid = (no_trans_sqrt-1)/2 * transducer_spacing;
    x_coord = linspace(xyz0(1) + x_mid - surf_side_mid, xyz0(1) + x_mid + surf_side_mid, no_trans_sqrt)';
    y_coord = linspace(xyz0(1) + y_mid - surf_side_mid, xyz0(1) + y_mid + surf_side_mid, no_trans_sqrt)';
    for x_c = 1:length(x_coord)
        for y_c = 1:length(y_coord)
            trans_pos(x_c + length(x_coord)*(y_c-1),:) = [x_coord(x_c) y_coord(y_c) d0];
        end
    end
    
end

%% Create forward operator A (generator)

padding = eir_length - 1;
H0_gen = Generate_Sparse_Delta_System_Matrix(Nx_gen-1, Ny_gen-1, 0, d_gen, Nt, dt, trans_pos, ...
    xyz0, [], [], [], 0);
eir_derivative = ConvolvedEIR(f0, FWHM, dt, eir_length)';
E = eir_derivative;

ff = @(x) D2D_Forward_Implicit(...
x, ...
H0_gen, ...
E, ...
no_trans, ...
padding);

A_gen = implicit2explicit(ff, (Nt+1)*no_trans, N_gen);

%% Create all speckled recordings

if ~strcmp(no_spec, 'inf')
    rec = zeros(size(A_gen,1), no_spec);
    rho = Rho(:);
    fprintf('Computing speckled recordings...\n')
    for s = 1:no_spec
        spec = Generate_2D_Speckle(Nx_gen, Ny_gen, r);
        rec(:, s) = A_gen*(rho.*spec(:));
        if mod(s,100) == 0
            fprintf('\tFinished computing speckled recording no. %d\n', s);
        end
    end

    % Add noise
    max_mag = max(abs(rec(:)));
    rec = rec + noise_level*max_mag*randn(size(rec));
end

%% Compute empirical estimates

if strcmp(no_spec, 'inf')
    R = diag(Rho(:));
    C_z_hat = A*R*C_E*R*A.' + C_eps;
else
    mu_z_hat = 1/no_spec * sum(rec,2);
    C_z_hat = rec*rec.';
    C_z_hat = C_z_hat/no_spec - mu_z_hat*mu_z_hat.';
end

%% Create forward operator A (reconstruction)

padding = eir_length - 1;
H0_rec = Generate_Sparse_Delta_System_Matrix(Nx_rec-1, Ny_rec-1, 0, d_rec, Nt, dt, trans_pos, ...
    xyz0, [], [], [], 0);
eir_derivative = ConvolvedEIR(f0, FWHM, dt, eir_length)';
E = eir_derivative;

ff = @(x) D2D_Forward_Implicit(...
x, ...
H0_rec, ...
E, ...
no_trans, ...
padding);

A_rec = implicit2explicit(ff, (Nt+1)*no_trans, N_rec);
C_eps_rec = sparse(eye(size(A_rec,1)));

%% Speckled reconstruct

C_eps_gen = (noise_level*max_mag)^2*sparse(eye(size(A_gen,1)));

tic
lambda = 1e+31;
lambda_2 = lambda;
I = eye(N_rec);
C_E_sqrt = sqrtm(C_E);
C_y_hat = C_z_hat - C_eps_rec;

dA = decomposition([A_rec; sqrt(lambda)*I]);
temp1 = dA \ [C_y_hat; zeros(N_rec, size(C_y_hat,2))];
temp2 = (dA \ [temp1.'; zeros(N_rec)]).';
temp3 = C_E_sqrt*temp2*C_E_sqrt;
S = sqrtm(temp3);

dB = decomposition([C_E_sqrt; sqrt(lambda_2)*I]);
temp4 = dB \ [S; zeros(N_rec)];
Sout = (dB \ [temp4.'; zeros(N_rec, size(temp4, 1))]).';
out = real(diag(Sout));

toc

%% Uniform reconstruction

lambda_unif = 1e+30;
if strcmp(no_spec, 'inf')
    avg_signal = A_gen*rho;
else
    avg_signal = mean(rec,2);
end
out_unif = [A_rec; sqrt(lambda_unif)*eye(N_rec)] \ [avg_signal; zeros(N_rec, size(avg_signal, 2))];

%% Plot and save stuff

if plotting || saving
    Out = reshape(out,Nx_rec,Ny_rec);
    Out_unif = reshape(out_unif,Nx_rec,Ny_rec);
end

if plotting
    f1 = figure;
    figure(f1)
    subplot(2,2,1)
    imagesc(Rho)
    title('True density')

    subplot(2,2,2)
    imagesc(Out)
    title('Speckled reconstruction')

    subplot(2,2,3)
    imagesc(Out_unif)
    title('Uniform reconstruction')

    subplot(2,2,4)
    imagesc(spec)
    title('Speckle example')
end

if saving
    fname = "N_gen_" + num2str(Nx_gen) + "_N_rec_" + num2str(Nx_rec) + "_r_" + num2str(r) + "_" + string(transducer_setting(1:4));
    save(fname, 'Rho', 'Out', 'Out_unif', 'spec')
    save("avg_signal_"+transducer_setting, 'avg_signal', 'trans_pos')
end
