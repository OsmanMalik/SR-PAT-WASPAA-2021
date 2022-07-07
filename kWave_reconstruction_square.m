% Before running this, first run waspaa_experiments.m for the square
% transducer geometry. This is necessary since that script creates the file
% 'avg_signal_surface array.mat' which is loaded in this script. See the
% README file for further details.
%
% This script requires k-Wave which can be downloaded at
% http://www.k-wave.org/.

% Load data 
load('avg_signal_surface array')

% Some settings
Nx_rec = 119;  % No. points along x-axis (reconstruction) 
Ny_rec = 119;  % No. points along y-axis (reconstruction) 
Nz_rec = 51;
d_rec = 4e-6*40/(80);  % Spatial grid spacing [m]  Default 4e-6
Nt = 200;
dt = 1e-9;
f0 = 50e6;
FWHM = 25e6;

% Grid
kgrid = kWaveGrid(Nx_rec, d_rec, Ny_rec, d_rec, Nz_rec, d_rec);
kgrid.setTime(Nt+1, dt);

% Medium
medium.sound_speed = 1500; %[m/s]

% Source
source.p0 = 0;

% Sensors
sensor.mask = [trans_pos(:,1:2) - mean(trans_pos(:,1)), trans_pos(:,3)].';
signal_reshape = reshape(avg_signal, Nt+1, size(trans_pos, 1));
sensor.time_reversal_boundary_data = signal_reshape.';
sensor.frequency_response = [f0 100*(FWHM/f0)];

% Other input arguments
input_args = {'PMLInside', false, 'Smooth', false, 'PlotPML', false,'PlotSim', false};

% Do reconstruction
p0_recon = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
p0_recon = p0_recon(20:100, 20:100, 26);

% Save results
save('kWave_recon_square', 'p0_recon')
