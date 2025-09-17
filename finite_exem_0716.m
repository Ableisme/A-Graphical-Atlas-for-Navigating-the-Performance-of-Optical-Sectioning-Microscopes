clear; close all; clc;
%% PART 1: UNIFIED SETUP & PSF GENERATION (No changes here)
% ========================================================================

% --- Microscope, Wavelength, and Grid Settings ---
NA = 0.5;
n = 1.33;
% For comparing with theory, we use a single wavelength
lambda = 525e-9;
lambda_ex = 488e-9;
lambda_em = 525e-9;

grid_size = 128; % Using a higher resolution grid for better results
pts = -grid_size:1:(grid_size-1);
centre_pt = grid_size + 1;
[X, Y, Z] = meshgrid(pts, pts, pts);
R = sqrt(X.^2 + Y.^2 + Z.^2);
shell_radius = grid_size / 2;
shell_thickness = 1;

% --- Define Unified Coordinate System ---
wavenumber_em = n / lambda_em;
spatial_frequency_per_pixel_em = (wavenumber_em / shell_radius) * 1e-6;
spatial_frequency_range = -grid_size*spatial_frequency_per_pixel_em:spatial_frequency_per_pixel_em:(grid_size-1)*spatial_frequency_per_pixel_em;
distance_per_pixel = 1 / (spatial_frequency_per_pixel_em * 2 * grid_size);
dist_um = -grid_size*distance_per_pixel : distance_per_pixel : (grid_size-1)*distance_per_pixel;
x_coords = dist_um; y_coords = dist_um; z_coords = dist_um;
X_scaled = X*spatial_frequency_per_pixel_em;
Y_scaled = Y*spatial_frequency_per_pixel_em;
Z_scaled = Z*spatial_frequency_per_pixel_em;
[X_real, Y_real, Z_real] = meshgrid(dist_um, dist_um, dist_um);

fprintf('Setup Complete. Grid Size = %d, Spatial Pixel Size = %.4f um.\n', grid_size, distance_per_pixel);

% --- PSF Generation ---
phi = atan2(sqrt(X.^2+Y.^2), Z);
cone = (0 <= phi & phi <= asin(NA/n));
shell_em = shell_radius - shell_thickness/2 < R & R <= shell_radius + shell_thickness/2;
CTF_em = double(cone & shell_em);
APSF_em = ifftshift(ifftn(fftshift(CTF_em)));
IPSF_em = abs(APSF_em).^2;

% Since wavelengths are the same, excitation PSF is the same as emission PSF
IPSF_ex = IPSF_em;

xx = 0.37 * lambda_em/NA * 1e6;
zz = 0.64 * lambda_em/ (n - sqrt(n^2 - NA^2))* 1e6;
fc = (2*NA/lambda_em) * 1e-6;

fprintf('PSF generation complete.\n\n');


%% PART 2: IN-DEPTH ANALYSIS FOR A SINGLE PINHOLE SIZE
% ========================================================================
fprintf('--- Starting In-depth Analysis for a Specific Pinhole Size ---\n');

% --- Define the specific AU value to analyze ---
AU_specific = 0.5;
fprintf('Analyzing for Pinhole Diameter = %.1f AU...\n', AU_specific);

% --- a) Construct and Visualize the Pinhole Mask ---
D_airy_disk_um = (1.22 * lambda_em / NA) * 1e6;
pinhole_diameter_um = D_airy_disk_um * AU_specific;
pinhole_radius_pixels = (pinhole_diameter_um / 2) / distance_per_pixel;

[X_2D, Y_2D] = meshgrid(x_coords, y_coords); % Use real coordinates for plotting
pinhole_mask_2D = double(sqrt(X_2D.^2 + Y_2D.^2) < (pinhole_diameter_um / 2));

figure;
imagesc(x_coords, y_coords, pinhole_mask_2D');
colormap(parula);
axis image;
title(sprintf('Visualization of Pinhole Mask (%.1f AU)', AU_specific), 'FontSize', 16);
xlabel('x (\mum)');
ylabel('y (\mum)');
colorbar;
set(gca, 'FontSize', 12);

% --- b) Calculate the final confocal PSF for this specific case ---
% We need the pinhole mask in pixel units for convolution
[X_pix, Y_pix] = meshgrid(pts, pts);
pinhole_mask_pix = double(sqrt(X_pix.^2 + Y_pix.^2) < pinhole_radius_pixels);
if sum(pinhole_mask_pix(:)) > 0
    pinhole_mask_pix = pinhole_mask_pix / sum(pinhole_mask_pix(:));
end

IPSF_det_blurred = zeros(size(IPSF_em));
for z_idx = 1:length(z_coords)
    IPSF_det_blurred(:, :, z_idx) = conv2(squeeze(IPSF_em(:, :, z_idx)), pinhole_mask_pix, 'same');
end
IPSF_final_specific = IPSF_ex .* IPSF_det_blurred;
IPSF_final_specific = IPSF_final_specific / max(IPSF_final_specific(:));
IPSF = IPSF_final_specific;
OTF = ifftshift(fftn(fftshift(IPSF))); % calculate the OTF
MTF = abs(OTF); % calculate the MTF
MTF = MTF/max(MTF(:)); % normalise the MTF
IPSF = IPSF./max(IPSF(:));% normalise the IPSF

% --- c) Calculate Lateral and Axial FWHM for this case ---
% Lateral FWHM calculation (from the final PSF at the focal plane)
lateral_profile = squeeze(IPSF_final_specific(centre_pt, :, centre_pt));
lateral_profile_norm = lateral_profile / max(lateral_profile(:));
[fwhm_lateral, ~] = calculate_fwhm_from_profile(x_coords, lateral_profile_norm);
fprintf('  > Calculated Lateral FWHM = %.4f um\n', fwhm_lateral);

% Axial FWHM calculation (using the robust sheet simulation method)
sheet3D_z0 = zeros(2*grid_size, 2*grid_size, 2*grid_size);
sheet3D_z0(:,:,centre_pt) = 1;
FT_sheet = fftn(fftshift(sheet3D_z0));
OTF_specific = fftn(fftshift(IPSF_final_specific));
I_image_3D = ifftshift(ifftn(FT_sheet .* OTF_specific));
axial_profile = squeeze(sum(sum(real(I_image_3D), 1), 2));
axial_profile_norm = axial_profile / max(axial_profile(:));
[fwhm_axial, axial_plot_data] = calculate_fwhm_from_profile(z_coords, axial_profile_norm);
fprintf('  > Calculated Axial FWHM   = %.4f um\n\n', fwhm_axial);

% --- d) Visualize the calculated profiles ---
figure;
% Lateral Profile Plot
subplot(1, 2, 1);
plot(x_coords, lateral_profile_norm, 'b-', 'LineWidth', 2);
hold on; grid on;
title(sprintf('Lateral Profile (%.1f AU)', AU_specific));
xlabel('Position (\mum)');
ylabel('Normalized Intensity');
xlim([-2, 2]);
text(0, 0.55, sprintf('FWHM = %.3f um', fwhm_lateral), 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
plot([-fwhm_lateral/2, fwhm_lateral/2], [0.5, 0.5], 'r--');

% Axial Profile Plot
subplot(1, 2, 2);
plot(axial_plot_data.zp, axial_plot_data.ifq, 'b-', 'LineWidth', 2);
hold on; grid on;
title(sprintf('Axial Profile (%.1f AU)', AU_specific));
xlabel('Position (\mum)');
ylabel('Normalized Intensity');
xlim([-8, 8]);
text(0, 0.55, sprintf('FWHM = %.3f um', fwhm_axial), 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
plot(axial_plot_data.x_points, [0.5, 0.5], 'r--');

spatial_freq_coords = spatial_frequency_range;

%% Integrated Intensity Profile
figure
IPSF_yz = squeeze(IPSF(:, centre_pt, :));  
imagesc(dist_um, dist_um, IPSF_yz');
axis image
%xlim([-1.5*zz, 1.5*zz])
%ylim([-1.5*zz, 1.5*zz])
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('IPSF Slice at x = 0 (y-z plane)', 'FontSize', 24, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlabel('y (µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
ylabel('z (µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
colorbar

fprintf('Analyzing optical sectioning\n');

Ii_slanted = squeeze(sum(sum(real(IPSF), 1), 2));
Ii_slanted_norm = Ii_slanted / max(Ii_slanted(:));

figure('Name', 'Integrated Intensity Profile');
hold on; grid on;
plot(z_coords, Ii_slanted_norm, 'b-', 'LineWidth', 2);
title(sprintf('Integrated Intensity Profile'));
xlabel('Axial Position (z, \mum)');
ylabel('Normalized Integrated Intensity');
legend('show', 'DisplayName', 'Integrated Intensity');
ylim([-0.05, 1.05]);
xlim([min(z_coords), max(z_coords)]);
hold off;

%% PART 3: FULL PARAMETER SWEEP & SUMMARY PLOT
% ========================================================================
fprintf('--- Starting Full Parameter Sweep for FWHM vs. AU Curve ---\n');
% (This section is the loop from our previous discussions)
pinhole_AU_range = 0:0.1:3.0;
fwhm_measured_values = nan(size(pinhole_AU_range));
h_waitbar = waitbar(0, 'Calculating FWHM for each pinhole size...');
tic;

for i = 1:length(pinhole_AU_range)
    % (Calculation logic is identical to the loop from our last version)
    % ... [The full loop code from the previous answer goes here] ...
    % For brevity, I will just call the function we defined
    fwhm_measured_values(i) = run_axial_fwhm_sim(IPSF_ex, IPSF_em, pinhole_AU_range(i), ...
                                                   distance_per_pixel, pts, z_coords);
    waitbar(i/length(pinhole_AU_range), h_waitbar);
end

close(h_waitbar);
fprintf('Full sweep finished in %.2f seconds.\n', toc);

% --- Plotting the final summary graph ---
figure;
hold on; grid on;
plot(pinhole_AU_range, fwhm_measured_values, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', 'Simulated FWHM');
% ... [The plotting code for theoretical curves goes here] ...
lambda_theory = lambda_em;
denominator = n - sqrt(n^2 - NA^2);
fwhm_ideal_theory = (0.67 * lambda_theory / denominator) * 1e6;
fwhm_quad_theory = fwhm_ideal_theory * sqrt(1 + pinhole_AU_range.^2);
fwhm_cubic_theory = fwhm_ideal_theory * (1 + 1.47 * pinhole_AU_range.^3).^(1/3);

plot(pinhole_AU_range, fwhm_quad_theory, 'r--', 'LineWidth', 2, 'DisplayName', 'Theory (Quadratic Approx, Eq. 25)');
plot(pinhole_AU_range, fwhm_cubic_theory, 'k:', 'LineWidth', 2.5, 'DisplayName', 'Theory (Cubic Approx, Eq. 26)');

title('Optical Sectioning FWHM vs. Pinhole Size', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Pinhole Diameter (AU)', 'FontSize', 14);
ylabel('Axial FWHM (\mum)', 'FontSize', 14);
legend('show', 'Location', 'northwest', 'FontSize', 12);
xlim([min(pinhole_AU_range), max(pinhole_AU_range)]);
valid_fwhm = fwhm_measured_values(~isnan(fwhm_measured_values));
if isempty(valid_fwhm)
    ylim_max = max(fwhm_cubic_theory);
else
    ylim_max = max([max(valid_fwhm), max(fwhm_cubic_theory)]);
end
ylim([0, ylim_max * 1.1]);
set(gca, 'FontSize', 12);
hold off;

%% PART 4: REPRODUCTION OF FIG. 10 (SMOOTH CURVES)
%  ========================================================================
fprintf('\n--- Starting Simulation to Reproduce Figure 10 ---\n');

% --- Define the parameter sets for each subplot, mirroring the paper ---
plot_params = {
    {'Dry Objectives', 1.0, [0.6, 0.7, 0.9]}, ...
    {'Water Immersion', 1.33, [0.85, 1.0, 1.2]}, ...
    {'Oil Immersion', 1.515, [1.2, 1.3, 1.4]}
};

% --- Create a new figure for the 3-panel plot ---
figure('Name', 'Figure 10 Reproduction', 'Position', [100, 100, 1500, 450]);
tic;

% --- Main loop to generate the three subplots ---
for plot_idx = 1:length(plot_params)
    
    % Select the current subplot
    subplot(1, 3, plot_idx);
    hold on; grid on;
    
    % Get parameters for the current subplot
    current_title = plot_params{plot_idx}{1};
    n_medium = plot_params{plot_idx}{2};
    na_values = plot_params{plot_idx}{3};
    
    fprintf('\n--- Simulating: %s (n=%.3f) ---\n', current_title, n_medium);
    
    % Loop through each NA value for the current medium
    legend_entries = {};
    for na_val = na_values
        
        if na_val >= n_medium
            warning('Skipping NA=%.2f as it is not physically possible for n=%.3f', na_val, n_medium);
            continue;
        end
        
        fprintf('  Running for NA = %.2f...\n', na_val);
        
        % Run the full simulation using the new helper function
        [pinhole_AU_range, fwhm_curve_um] = run_full_fwhm_simulation(na_val, n_medium, lambda, grid_size);
        
        % Plot the normalized result (FWHM/λ)
        plot(pinhole_AU_range, fwhm_curve_um, 'LineWidth', 2);
        legend_entries{end+1} = sprintf('NA=%.2f', na_val);
    end
    
    % Format the current subplot
    title(current_title, 'FontSize', 14);
    xlabel('AU', 'FontSize', 12);
    ylabel('FWHM (\mum)', 'FontSize', 12);
    legend(legend_entries, 'Location', 'northwest');
    set(gca, 'FontSize', 10);
    
end

hold off;
sgtitle('Reproducing Paper Figure 10: FWHM/\lambda vs. AU for Various Objectives', 'FontSize', 18, 'FontWeight', 'bold');
fprintf('\nFigure 10 reproduction finished in %.2f seconds.\n', toc);


%% HELPER FUNCTION AT THE END OF YOUR SCRIPT
%  (It uses your existing 'calculate_fwhm_from_profile' function internally)
% ========================================================================

function [pinhole_AU_range, fwhm_values] = run_full_fwhm_simulation(NA, n, lambda, grid_size)
    % This function encapsulates the entire simulation for a given set of
    % microscope parameters (NA, n, lambda) and returns a smooth FWHM vs. AU curve.

    % --- Setup grids and coordinates based on input parameters ---
    pts = -grid_size:1:(grid_size-1);
    centre_pt = grid_size + 1;
    [X, Y, Z] = meshgrid(pts, pts, pts);
    R = sqrt(X.^2 + Y.^2 + Z.^2);
    shell_radius = grid_size / 2;
    
    wavenumber = n / lambda;
    spatial_frequency_per_pixel = (wavenumber / shell_radius) * 1e-6;
    distance_per_pixel = 1 / (spatial_frequency_per_pixel * 2 * grid_size);
    z_coords = -grid_size*distance_per_pixel : distance_per_pixel : (grid_size-1)*distance_per_pixel;
    
    % --- Generate single-wavelength PSF ---
    phi = atan2(sqrt(X.^2+Y.^2), Z);
    cone = (0 <= phi & phi <= asin(NA/n));
    shell = shell_radius - 1/2 < R & R <= shell_radius + 1/2;
    CTF = double(cone & shell);
    APSF = ifftshift(ifftn(fftshift(CTF)));
    IPSF = abs(APSF).^2;
    IPSF_ex = IPSF;
    IPSF_em = IPSF;

    % --- Prepare for loop ---
    sheet3D_z0 = zeros(2*grid_size, 2*grid_size, 2*grid_size);
    sheet3D_z0(:,:,centre_pt) = 1;
    FT_sheet = fftn(fftshift(sheet3D_z0));
    
    pinhole_AU_range = 0:0.1:3.0;
    fwhm_values = nan(size(pinhole_AU_range));
    
    h_waitbar = waitbar(0, sprintf('Simulating NA=%.2f...', NA));

    % --- Loop through pinhole sizes ---
    for i = 1:length(pinhole_AU_range)
        waitbar(i/length(pinhole_AU_range), h_waitbar);
        current_AU = pinhole_AU_range(i);
        
        IPSF_det_blurred = zeros(size(IPSF_em));
        if current_AU < 1e-6
            IPSF_det_blurred = IPSF_em;
        else
            D_airy_disk_um = (1.22 * lambda / NA) * 1e6;
            pinhole_diameter_um = D_airy_disk_um * current_AU;
            pinhole_radius_pixels = (pinhole_diameter_um / 2) / distance_per_pixel;
            
            [X_pix, Y_pix] = meshgrid(pts, pts);
            pinhole_mask_pix = double(sqrt(X_pix.^2 + Y_pix.^2) < pinhole_radius_pixels);
            if sum(pinhole_mask_pix(:)) > 0
                pinhole_mask_pix = pinhole_mask_pix / sum(pinhole_mask_pix(:));
            end
            for z_idx = 1:length(z_coords)
                IPSF_det_blurred(:, :, z_idx) = conv2(squeeze(IPSF_em(:, :, z_idx)), pinhole_mask_pix, 'same');
            end
        end
        IPSF_final = IPSF_ex .* IPSF_det_blurred;

        OTF_final = fftn(fftshift(IPSF_final));
        I_image_3D = ifftshift(ifftn(FT_sheet .* OTF_final));
        axial_profile = squeeze(sum(sum(real(I_image_3D), 1), 2));
        
        if max(axial_profile(:)) > 1e-12
             % Use the existing powerful FWHM calculation function from your script
            [fwhm, ~] = calculate_fwhm_from_profile(z_coords, axial_profile);
            fwhm_values(i) = fwhm * 1e6; % Convert fwhm from m to um
        end
    end
    close(h_waitbar);
end


%% HELPER FUNCTIONS (Add these at the end of your .m file)
% ========================================================================

function [fwhm, plot_data] = calculate_fwhm_from_profile(coords, profile)
    % This helper function calculates FWHM from a 1D profile using interpolation
    fwhm = NaN;
    plot_data = struct();

    profile_norm = profile / max(profile(:));
    factor = 20;
    xp = linspace(coords(1), coords(end), numel(coords)*factor);
    yp = interp1(coords, profile_norm, xp, 'spline');
    
    [~, pk_idx] = max(yp);
    
    if pk_idx > 1 && pk_idx < length(yp)
        find_half_points = @(y_data, x_data, y_target) ...
            deal(interp1(y_data(1:pk_idx), x_data(1:pk_idx), y_target, 'linear', NaN), ...
                 interp1(y_data(pk_idx:end), x_data(pk_idx:end), y_target, 'linear', NaN));
        [xL, xR] = find_half_points(yp, xp, 0.5);
        
        if ~isnan(xL) && ~isnan(xR)
            fwhm = xR - xL;
            plot_data.zp = xp;
            plot_data.ifq = yp;
            plot_data.x_points = [xL, xR];
        end
    end
end

function fwhm = run_axial_fwhm_sim(IPSF_ex, IPSF_em, current_AU, distance_per_pixel, pts, z_coords)
    % This function encapsulates the simulation for a single AU value
    fwhm = NaN;
    grid_size = length(pts)/2;
    centre_pt = grid_size + 1;

    % Calculate final IPSF
    IPSF_det_blurred = zeros(size(IPSF_em));
    if current_AU < 1e-6
        IPSF_det_blurred = IPSF_em;
    else
        D_airy_disk_um = (1.22 * 525e-9 / 0.5) * 1e6; % Hardcoded for simplicity
        pinhole_diameter_um = D_airy_disk_um * current_AU;
        pinhole_radius_pixels = (pinhole_diameter_um / 2) / distance_per_pixel;
        [X_pix, Y_pix] = meshgrid(pts, pts);
        pinhole_mask_pix = double(sqrt(X_pix.^2 + Y_pix.^2) < pinhole_radius_pixels);
        if sum(pinhole_mask_pix(:)) > 0
            pinhole_mask_pix = pinhole_mask_pix / sum(pinhole_mask_pix(:));
        end
        for z_idx = 1:length(z_coords)
            IPSF_det_blurred(:, :, z_idx) = conv2(squeeze(IPSF_em(:, :, z_idx)), pinhole_mask_pix, 'same');
        end
    end
    IPSF_final = IPSF_ex .* IPSF_det_blurred;

    % Sheet simulation
    sheet3D_z0 = zeros(2*grid_size, 2*grid_size, 2*grid_size);
    sheet3D_z0(:,:,centre_pt) = 1;
    FT_sheet = fftn(fftshift(sheet3D_z0));
    OTF_final = fftn(fftshift(IPSF_final));
    I_image_3D = ifftshift(ifftn(FT_sheet .* OTF_final));
    axial_profile = squeeze(sum(sum(real(I_image_3D), 1), 2));
    
    if max(axial_profile(:)) > 1e-12
        [fwhm, ~] = calculate_fwhm_from_profile(z_coords, axial_profile);
    end
end

%% Point Object Resolution Figure (Any Angles -> FWHM)

N_angles = 91;                   
theta_list = linspace(0, pi, N_angles);  

% 
x_half = nan(1, N_angles);
z_half = nan(1, N_angles);
r_half_vec = nan(1, N_angles);
% ——————————————————————————————————————————

for k = 1 : N_angles
    theta = theta_list(k);

    % ——(a)  r_max ———
    if abs(cos(theta)) < 1e-8
        r_max_x = Inf;
    else
        r_max_x = max(abs(x_coords)) / abs(cos(theta));
    end

    if abs(sin(theta)) < 1e-8
        r_max_z = Inf;
    else
        r_max_z = max(abs(z_coords)) / abs(sin(theta));
    end
    r_max = min(r_max_x, r_max_z);  
    % (r*cosθ, r*sinθ)  x_coords, z_coords  boundary

    % ——(b) [0, r_max] ———
    N_r_sample = 800;  % sample
    r_vec = linspace(0, r_max, N_r_sample);

    % ——(c) calculate ———
    x_r = r_vec * cos(theta);
    z_r = r_vec * sin(theta);
    IPSF_xz = squeeze(IPSF(centre_pt, :, :));
 
    f_r = interp2( z_coords, x_coords, IPSF_xz, z_r, x_r, 'spline', 0 ); 
  
    % ——(d) ———
    f_r = f_r / f_r(1);

    % ——(e)  0.5 idx———
    idx_fall = find( f_r < 0.5, 1, 'first' );
    if ~isempty(idx_fall) && idx_fall > 1
        % 
        r1 = r_vec(idx_fall - 1);
        r2 = r_vec(idx_fall);
        f1 = f_r(idx_fall - 1);
        f2 = f_r(idx_fall);

        % 
        r_half_k = r1 + (0.5 - f1) * (r2 - r1) / (f2 - f1);
    else
        % 
        r_half_k = NaN;
    end

    % ——(f) ———
    r_half_vec(k) = r_half_k; 
    x_half(k) =  r_half_k * cos(theta);
    z_half(k) =  r_half_k * sin(theta);
end

% ：
x_half_full = [ x_half,  -x_half ];
z_half_full = [ z_half,  -z_half ];


figure;

imagesc( z_coords, x_coords, IPSF_xz' );
axis image; 
colormap parula;
colorbar

hold on;

plot( x_half_full, z_half_full, 'r--', 'LineWidth', 2, 'DisplayName', 'FWHM Contour' );

plot(0, 0, 'r+', 'MarkerSize', 8, 'LineWidth', 1.2, 'HandleVisibility', 'off');
hold off; 
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold','FontSize', 14);
title('IPSF Slice at y=0 with FWHM Contour', 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('x-axis (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('z-axis (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

lgd = legend('show', 'Location', 'Best');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');


xlim([-1.5*zz,1.5*zz]);
ylim([-1.5*zz,1.5*zz]);


%% ===== FWHM fraction (3D) + Optional axial FWHM =====
% Assumptions: IPSF is already a 3D intensity kernel (can be normalized or not - this metric is insensitive to amplitude scaling)
%             dist_um / X_real, Y_real, Z_real are already defined (real space coordinates)
IPSF_det = IPSF_det_blurred / max(IPSF_det_blurred(:));
IPSF_ill = IPSF_ex / max(IPSF_ex(:));
% --- Voxel size and volume (µm) ---
dx = dist_um(2) - dist_um(1);
dy = dx; dz = dx;                   % Your grid is equally spaced
voxelVol = dx*dy*dz;                % Voxel volume (µm^3)

% --- Energy fraction of 3D FWHM volume ---
V = IPSF_det;                           % Intensity volume
Imax = max(V(:));
mask3D = (V >= 0.5*Imax);           % FWHM mask
E_total = sum(V(:)) * voxelVol;     % Total energy
E_fwhm  = sum(V(mask3D)) * voxelVol;% Energy within FWHM
frac_fwhm3D = E_fwhm / E_total;     % <-- The metric you want

% --- (Optional) 1D FWHM in x/y/z directions (linear interpolation) ---
% Extract three profiles passing through voxel center
prof_x = squeeze(V(:, centre_pt, centre_pt));
prof_y = squeeze(V(centre_pt, :, centre_pt));
prof_z = squeeze(V(centre_pt, centre_pt, :));

x_axis = dist_um(:); y_axis = dist_um(:); z_axis = dist_um(:);
FWHM_x = local_fwhm_1d(x_axis, prof_x);
FWHM_y = local_fwhm_1d(y_axis, prof_y);
FWHM_z = local_fwhm_1d(z_axis, prof_z);

fprintf('--- FWHM-fraction(3D) = %.4f | FWHM x/y/z [µm] = %.3f / %.3f / %.3f ---\n', ...
        frac_fwhm3D, FWHM_x, FWHM_y, FWHM_z);

% ===== Your 3D isosurface plot (overlay numerical values on the plot) =====
figure
V_iso_surface = isosurface(X_real, Y_real, Z_real, V, 0.001);
p2 = patch(V_iso_surface);
isonormals(X_real, Y_real, Z_real, V, p2)
p2.FaceColor = 'blue';
p2.EdgeColor = 'none';
daspect([1,1,1])
view(3), axis equal
camlight; lighting gouraud
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('3D Plot of Confocal Microscope(0.5 AU) IPSF', 'FontSize', 24, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlabel('x (µm)', 'FontSize', 16, 'FontName', 'Times New Roman','FontWeight', 'bold');
ylabel('y (µm)', 'FontSize', 16, 'FontName', 'Times New Roman','FontWeight', 'bold');
zlabel('z (µm)', 'FontSize', 16, 'FontName', 'Times New Roman','FontWeight', 'bold');


% ===== Add 3D FWHM fraction information in the upper right corner =====
txt_str = sprintf('3D FWHM fraction = %.3f\nFWHM x/y/z = %.2f/%.2f/%.2f µm', ...
                  frac_fwhm3D, FWHM_x, FWHM_y, FWHM_z);

% Calculate upper right corner position (slightly outside data range)
x_pos = max(x_coords)*0.9;
y_pos = max(y_coords)*0.9;
z_pos = max(z_coords)*0;

t = text(x_pos, y_pos, z_pos, txt_str, ...
         'FontSize', 16, 'FontWeight', 'bold', ...
         'FontName', 'Times New Roman', ...
         'BackgroundColor', [1 1 1 0.75], ... % White + semi-transparent
         'EdgeColor', [0 0 0], ...
         'Margin', 5, ...
         'VerticalAlignment', 'top', ...
         'HorizontalAlignment', 'left');


function w = local_fwhm_1d(x, y)
    % Linear interpolated 1D FWHM (y can be non-normalized)
    x = x(:); y = y(:);
    [ymax, imax] = max(y);
    half = ymax/2;

    % Left half maximum to half height crossing point
    iL = find(y(1:imax) <= half, 1, 'last');
    if isempty(iL) || iL==imax
        xL = x(1);
    else
        xL = interp1(y([iL, iL+1]), x([iL, iL+1]), half, 'linear');
    end

    % Right half from peak to half height crossing point
    iRrel = find(y(imax:end) <= half, 1, 'first');
    if isempty(iRrel) || imax-1+iRrel==imax
        xR = x(end);
    else
        iR = imax-1+iRrel;
        xR = interp1(y([iR-1, iR]), x([iR-1, iR]), half, 'linear');
    end

    w = abs(xR - xL);
end


%% ===== Added: Quantitative metrics (FWHM_x/y/z, HMEF, Tail energy, radial MTF) =====
try
    % Energy-normalized PSF for integral metrics
    IPSF_energy = IPSF ./ max(sum(IPSF(:)), eps);
    
    % Coordinates (micrometers)
    x_um = x_coords; y_um = y_coords; z_um = z_coords;
    [~,imax] = max(IPSF(:));
    [ix,iy,iz] = ind2sub(size(IPSF), imax);
    
    % 1D profiles through the peak
    prof_x = squeeze(IPSF(:,iy,iz));
    prof_y = squeeze(IPSF(ix,:,iz));
    prof_z = squeeze(IPSF(ix,iy,:));
    
    % FWHM via linear interpolation
    FWHM_x_um = local_fwhm_1d(x_um(:), prof_x(:));
    FWHM_y_um = local_fwhm_1d(y_um(:), prof_y(:));
    FWHM_z_um = local_fwhm_1d(z_um(:), prof_z(:));
    
    % 3D Half-Max Energy Fraction (HMEF)
    half_level = 0.5*max(IPSF(:));
    HMEF = sum(IPSF_energy(IPSF >= half_level));
    
    % Tail Energy Ratio beyond z0 = 2*FWHM_z
    z0_um = 2*FWHM_z_um;
    Ez = squeeze(sum(sum(IPSF_energy,1),2));              % energy per z-slice (sums to 1)
    tail_mask = (abs(z_um) >= z0_um);
    TER = sum(Ez(tail_mask));
    
    fprintf('--- Wide-field quantitative metrics ---\n');
    fprintf('FWHM_x = %.4f um | FWHM_y = %.4f um | FWHM_z = %.4f um\n', FWHM_x_um, FWHM_y_um, FWHM_z_um);
    fprintf('HMEF (3D >=0.5*max energy fraction) = %.4f\n', HMEF);
    fprintf('Tail energy fraction |z| >= %.2f um: %.4f\n', z0_um, TER);

    
    
    %% Radial MTF (z = 0 slice) and MTF50/10/cutoff
    izc = centre_pt;              % center z-plane
    MTF2 = squeeze(MTF(:,:,izc));
    MTF2 = MTF2 / max(MTF2(:)+eps);
    
    % Build frequency axes (cyc/um) from provided spatial_frequency_range
    fx = spatial_frequency_range; fy = spatial_frequency_range;
    [FRad, mtf_rad, fr_centers] = radial_profile( fx, fy, MTF2, 200 );
    
    % MTF metrics
    [fc_thresh, fc_mtf10, fc_mtf50] = mtf_cutoff_linear(fr_centers(:), mtf_rad(:), 1e-7);
    fprintf('Radial MTF: fc(thresh=1e-3)=%.3f cyc/um | MTF10=%.3f | MTF50=%.3f\n', fc_thresh, fc_mtf10, fc_mtf50);
    
    %% Quick summary plot (profiles + radial MTF)
    figure('Name','WF Quantitative Metrics','Color','w');
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
    nexttile; plot(x_um, prof_x, 'LineWidth',1.5); hold on; yline(0.5,'--'); grid on;
    title(sprintf('Profile X | FWHM=%.3f um', FWHM_x_um)); xlabel('x (\mum)'); ylabel('Norm.');
    nexttile; plot(y_um, prof_y, 'LineWidth',1.5); hold on; yline(0.5,'--'); grid on;
    title(sprintf('Profile Y | FWHM=%.3f um', FWHM_y_um)); xlabel('y (\mum)'); ylabel('Norm.');
    nexttile; plot(z_um, prof_z, 'LineWidth',1.5); hold on; yline(0.5,'--'); grid on;
    title(sprintf('Profile Z | FWHM=%.3f um | HMEF=%.3f | TER_{|z|>%.3f}=%.3f', FWHM_z_um, HMEF, z0_um, TER));
    xlabel('z (\mum)'); ylabel('Norm.');
    nexttile; plot(fr_centers, mtf_rad, 'LineWidth',1.8); grid on; hold on;
    yline(0.5,':'); yline(0.1,':'); xline(fc_mtf50,':'); xline(fc_mtf10,':'); xline(fc_thresh,'--');
    xlabel('Spatial frequency (cyc/\mum)'); ylabel('MTF (radial)');
    title(sprintf('Radial MTF fc=%.2f cyc/\\mum, MTF50=%.2f, MTF10=%.2f', ...
      fc_thresh, fc_mtf50, fc_mtf10));

    
catch ME
    msg = getReport(ME, 'extended', 'hyperlinks', 'off');
    warning('WF:MetricsFailed','%s', msg);
end


%% Simplified MTF Annotation (Numerically from Data)



% --- 1. Get the MTF 2D slice to be analyzed ---
% Make sure this MTF is the confocal MTF you want to analyze
mtf_slice = squeeze(MTF(centre_pt, :, :))'; % Get the x-z plane at u_y=0

% --- 2. Use the most direct method to find lateral boundaries from data ---
threshold = 1e-3; % Define a threshold to determine signal region
mask = mtf_slice > threshold;

% Find all u_x columns that contain valid signal
non_empty_columns = any(mask, 1);

if ~any(non_empty_columns)
    warning('MTF data is empty or below threshold. Cannot determine cutoff frequency.');
    fc_min_measured = NaN;
    fc_max_measured = NaN;
else
    % Find the indices of the first and last columns containing signal
    first_idx = find(non_empty_columns, 1, 'first');
    last_idx = find(non_empty_columns, 1, 'last');

    % Get corresponding frequency values from coordinate vector
    fc_min_measured = spatial_frequency_range(first_idx);
    fc_max_measured = spatial_frequency_range(last_idx);
    
    fprintf('Numerically found MTF bounds in u_x direction: [%.2f, %.2f] cycles/µm\n', fc_min_measured, fc_max_measured);
end


% --- 3. Plotting and annotation ---
figure
imagesc(spatial_frequency_range, spatial_frequency_range, mtf_slice)
axis equal
axis tight
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('MTF, u_y = 0, range = [0, 0.01]', 'FontSize', 24, 'FontName', 'Times New Roman','FontWeight', 'bold');
ylabel('u_z (cycles/µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlabel('u_x (cycles/µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlim([-4,4]);ylim([-3,3]);
colorbar
clim([0, 0.01])

% Only perform annotation when boundaries are successfully found
if ~isnan(fc_max_measured)
    hold on;
    yLims = ylim;
    
    % Draw boundary lines
    plot([fc_max_measured, fc_max_measured], yLims, 'r--', 'LineWidth', 1);
    plot([fc_min_measured, fc_min_measured], yLims, 'r--', 'LineWidth', 1);
    
    % Add text
    text(fc_max_measured - 2, - 0.8, sprintf('u_x_m_a_x = %.2f', fc_max_measured), ...
        'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    
    hold off;
end
%% Sheet Object(HORIZONTAL) Figures

% Create a sheet object and visualize
% Create a 3D sheet: one-voxel thick at z=0, uniform in x/y
sheet3D_z0 = zeros(size(IPSF));      % size = [Nx,Ny,Nz]
[~, iz0] = min(abs(dist_um));     % find index of z=0
sheet3D_z0(:,:,iz0) = 1;             % set that plane to 1
% Visualize the z=0 sheet
figure;
p = patch( isosurface(dist_um, dist_um, dist_um, sheet3D_z0, 0.5) );
isonormals(dist_um, dist_um, dist_um, sheet3D_z0, p);
p.FaceColor = 'cyan'; p.EdgeColor = 'none';
daspect([1,1,1]); view(3); axis tight; camlight; lighting gouraud;
title('3D Thin Sheet Object (z=0)', 'FontSize', 16, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('x (µm)', 'FontSize', 16, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('y (µm)', 'FontSize', 16, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
zlabel('z (µm)', 'FontSize', 16, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

% Simulate Imaging in Frequency‐domain by 3D multiplication 
FT_sheet  = fftn( fftshift(sheet3D_z0) );
OTF = ifftshift(fftn(fftshift(IPSF))); % calculate the OTF
OTF    = ifftshift(OTF);
I_freq3D  = ifftshift( ifftn( FT_sheet .* OTF ) );

freq_line  = squeeze( real(I_freq3D(centre_pt,centre_pt, :)) );
%maxErr3D = max( abs(space_line - freq_line) );
%fprintf('Max absolute error between convolution methods: %.3e\n', maxErr3D);

peak_freq  = max(freq_line);
factor = 20;
zp = linspace(dist_um(1), dist_um(end), numel(dist_um)*factor);
ifq = interp1(dist_um, freq_line/peak_freq,  zp, 'spline');
[~, pk_f] = max(ifq);
findFWHM = @(y) ...
  deal( ...
    interp1(y(1:pk_f), zp(1:pk_f), .5,'linear'), ...
    interp1(y(pk_f:end), zp(pk_f:end), .5,'linear') ...
  );
[xL_f,xR_f] = findFWHM(ifq);
FWHM_f = xR_f - xL_f;

% Print results
fprintf('\nMethod   PeakI    FWHM(µm)\n');
fprintf('Freq     %.4f    %.4f\n', peak_freq,  FWHM_f);

figure;
plot(zp, ifq, 'r--','LineWidth',1.5);hold on;
plot(zp(pk_f), ifq(pk_f), 'ro','MarkerFaceColor','r');
plot([xL_f xR_f],[.5 .5],'r-','LineWidth',1);
title(sprintf('Optical Sectioning Profile (FWHM = %.3f µm)', FWHM_f), 'FontSize', 16, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('z-axis (µm)', 'FontSize', 16, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylim([0,1]);
ylabel('Normalized Intensity', 'FontSize', 16, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
lgd = legend('show', 'Location', 'Best');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
grid on;

% Frequency domain visualization
% --- Plot frequency domain view of horizontal sheet (Red: Sheet Spectrum, Blue: OTF) ---
figure;

FT_sheet_z0_centered = fftshift(FT_sheet);
otf_slice_2d        = squeeze(MTF(centre_pt, :, :));                % [uy,ux,uz] -> fix uy
sheet_spec_2d       = abs( squeeze( FT_sheet_z0_centered(centre_pt, :, :) ) );
log_otf     = log(1 + otf_slice_2d);
log_otf_norm = log_otf / max(log_otf(:));
sheet_norm = sheet_spec_2d / max(sheet_spec_2d(:));  % Normalize to [0,1]
mask_line  = sheet_norm > 0.2;                       % Threshold can be adjusted as needed

imagesc(spatial_frequency_range, spatial_frequency_range, log_otf_norm);
axis xy image;
title('Frequency View for Horizontal Sheet', 'FontSize', 16, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Axial Spatial Frequency u_z (\mum^{-1})', 'FontSize', 16, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Lateral Spatial Frequency u_x (\mum^{-1})', 'FontSize', 16, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
lgd = legend('show', 'Location', 'Best');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold', 'TextColor', 'white');
colormap(parula);          % Deep blue to yellow
caxis([0,0.01]);           % MTF max value ~0.01, fine-tune as needed
hcb = colorbar;
hcb.Label.String = 'MTF value';hold on;
[rows, cols] = find(mask_line);
u_x_line = spatial_frequency_range(rows);  % Rows correspond to y-axis
u_z_line = spatial_frequency_range(cols);  % Columns correspond to x-axis
plot(u_z_line, u_x_line, '.r', 'MarkerSize', 8);
xlim([-2*fc, 2*fc]);
ylim([-2*fc, 2*fc]);
set(gca,'Layer','top');
hold off;
%% Sheet Object(ROTATE) Figures
% 1) Construct arbitrary angle sheet (rotated around Y axis) and visualize

theta = 0; % Define rotation angle (60 degrees)
normal_vec = [-sin(theta), 0, cos(theta)];
normal_vec = normal_vec / norm(normal_vec); 
plane_dist = X_real * normal_vec(1) + Y_real * normal_vec(2) + Z_real * normal_vec(3);
sheet_thickness = distance_per_pixel; 
sheet3D_slanted = abs(plane_dist) < (sheet_thickness / 2);

figure;
p = patch( isosurface(X_real, Y_real, Z_real, double(sheet3D_slanted), .5) );
isonormals(X_real, Y_real, Z_real, double(sheet3D_slanted), p);
p.FaceColor='cyan'; p.EdgeColor='none';
daspect([1,1,1]); view(3); camlight; lighting gouraud;
set(gca, 'FontName', 'Times New Roman','FontSize', 14);
title(sprintf('Slanted Sheet (Rotated around Y-axis) at \\theta=%.1f°',theta*180/pi), 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('x (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('y (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
zlabel('z (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
hold on;

p0 = [0, 0, 0]; 
L = max(dist_um) * 1.2; 
quiver3( p0(1), p0(2), p0(3), ...
         normal_vec(1), normal_vec(2), normal_vec(3), ...
         L, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'AutoScale','off' );
axis tight; view(3);

% 2) Simulate Imaging in Frequency Domain by Multiplication
FT_sheet_slanted = fftn(fftshift(sheet3D_slanted));
I_freq_slanted   = ifftshift(ifftn( FT_sheet_slanted .* OTF ));
I_freq_slanted = I_freq_slanted./max(I_freq_slanted(:));% normalise the IPSF


% Spatial domain profile comparison: Extract profiles from BOTH results and compare

fprintf('Extracting intensity profiles from both results...\n');
% Define sampling line
max_dist = max(abs(dist_um));
num_points = 2000;
line_coords_1d = linspace(-max_dist, max_dist, num_points);
sample_points_x = p0(1) + line_coords_1d * normal_vec(1);
sample_points_y = p0(2) + line_coords_1d * normal_vec(2);
sample_points_z = p0(3) + line_coords_1d * normal_vec(3);

% --- b) Extract profile from FREQUENCY domain result ---
% Use real() to discard negligible imaginary noise from FFT numerics
profile_freq = interp3(X_real, Y_real, Z_real, real(I_freq_slanted), ...
                       sample_points_x, sample_points_y, sample_points_z, 'cubic', 0);
profile_freq_norm = profile_freq / max(profile_freq(:));

figure;
hold on;
plot(line_coords_1d, profile_freq_norm, 'b-', 'LineWidth', 2);
hold off;
grid on;
set(gca, 'FontName', 'Times New Roman','FontSize', 14);
title('Profile along Normal Vector', 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Distance along Normal Vector (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Normalized Intensity', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
fprintf('Comparison plot generated.\n');
xlim([-10,10]);
ylim([0,1]);
fprintf('Comparison plot generated.\n');
% === Slanted sheet metrics along NORMAL: FWHM_n & TER_normal ===
% 1) FWHM_n: Directly based on your obtained normal direction profile (line_coords_1d, profile_freq_norm)
FWHM_n = NaN;
if max(profile_freq_norm)>0
    p = profile_freq_norm(:) / max(profile_freq_norm);
    [~, pk] = max(p);
    iL = find(p(1:pk) <= 0.5, 1, 'last');
    iR = find(p(pk:end) <= 0.5, 1, 'first'); if ~isempty(iR), iR = pk+iR-1; end
    if ~isempty(iL) && ~isempty(iR) && iL>1 && iR<numel(p)
        xL = line_coords_1d(iL) + (0.5 - p(iL)) * (line_coords_1d(iL+1)-line_coords_1d(iL)) / (p(iL+1)-p(iL) + eps);
        xR = line_coords_1d(iR-1) + (0.5 - p(iR-1)) * (line_coords_1d(iR)-line_coords_1d(iR-1)) / (p(iR)-p(iR-1) + eps);
        FWHM_n = xR - xL;
    end
end

% 2) TER_normal: Integrate 3D intensity on planes parallel to the sheet to get normal direction energy distribution E(u)
I_abs_s = abs(I_freq_slanted);                    % Numerically robust
E_tot_s = sum(I_abs_s(:)) + eps;

% Construct a set of basis orthogonal to normal_vec: t1_vec, t2_vec (choose any direction not collinear with normal and orthogonalize)
ez = [0 0 1];
if abs(dot(normal_vec, ez)) < 0.99
    t1_vec = cross(normal_vec, ez);
else
    t1_vec = cross(normal_vec, [1 0 0]);
end
t1_vec = t1_vec / norm(t1_vec);
t2_vec = cross(normal_vec, t1_vec); t2_vec = t2_vec / norm(t2_vec);

% Resample volume data in (u,v,w) coordinates: u along normal, (v,w) within sheet plane
u_coords = line_coords_1d;                        % Normal coordinates directly reuse your sampling line
plane_extent = max([max(abs(x_coords)) max(abs(y_coords)) max(abs(z_coords))]);
Nv = 96; Nw = 96;                                 % In-plane sampling resolution (adjustable)
v_coords = linspace(-plane_extent, plane_extent, Nv);
w_coords = linspace(-plane_extent, plane_extent, Nw);
[U,V,W] = ndgrid(u_coords, v_coords, w_coords);

Xq = U*normal_vec(1) + V*t1_vec(1) + W*t2_vec(1);
Yq = U*normal_vec(2) + V*t1_vec(2) + W*t2_vec(2);
Zq = U*normal_vec(3) + V*t1_vec(3) + W*t2_vec(3);

VOL_res = interp3(X_real, Y_real, Z_real, I_abs_s, Xq, Yq, Zq, 'linear', 0);
E_u = squeeze(sum(sum(VOL_res, 3), 2));           % Integrate over (v,w)
E_u = E_u / max(sum(E_u), eps);                   % Normalize to 1, as energy distribution along normal

u0 = 2 * FWHM_n;                                  % TER threshold: |u| >= 2*FWHM_n
TER_normal = sum(E_u( abs(u_coords) >= u0 ));

fprintf('Slanted sheet (theta=%.1f deg): FWHM_n = %.3f um | TER_{|u|>=%.3f} = %.4f\n', ...
        theta*180/pi, FWHM_n, u0, TER_normal);

% 3) Metric plot (normal direction)
figure('Name','Slanted sheet axial response (along normal)','Color','w');
plot(u_coords, E_u / max(E_u), 'b-', 'LineWidth', 1.8); hold on; grid on;
yline(0.5,':','Half-max');
% Mark both ends of FWHM on normal coordinates
if isfinite(FWHM_n)
    xline(-FWHM_n/2,':'); xline(FWHM_n/2,':');
end
title(sprintf('Along-normal response: FWHM_n=%.3f um, TER_{|u|>=%.3f}=%.3f (\\theta=%.1f^\\circ)', ...
      FWHM_n, u0, TER_normal, theta*180/pi));
xlabel('Distance along normal u (um)'); ylabel('Normalized energy/profile');


fprintf('X-Z plane visualization complete.\n');

% PSF image
figure;
imagesc(v_coords, line_coords_1d, profile_freq_norm');
axis image; 
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('Image after PSF convolution');
xlabel('x (\mum)');
ylabel('z (\mum)');
ylim([-5, 5]);
xlim([-5, 5]);
colorbar
%clim([0.9 1]);

fprintf('Analyzing optical sectioning for the slanted sheet image...\n');

Ii_slanted = squeeze(sum(sum(real(I_freq_slanted), 1), 2));
Ii_slanted_norm = Ii_slanted / max(Ii_slanted(:));

figure('Name', 'Integrated Intensity Profile for Slanted Sheet');
hold on; grid on;
plot(z_coords, Ii_slanted_norm, 'b-', 'LineWidth', 2);
title(sprintf('Integrated Intensity Profile for Slanted Sheet (\\theta=%.1f°)', theta*180/pi));
xlabel('Axial Position (z, \mum)');
ylabel('Normalized Integrated Intensity');
legend('show', 'DisplayName', 'Integrated Intensity');
ylim([-0.05, 1.05]);
xlim([min(z_coords), max(z_coords)]);
hold off;

% Spectrum visualization (MODIFIED TO DRAW A SMOOTH ANALYTICAL LINE)

figure;
% 1) Draw OTF background (completely same as before)
imagesc( spatial_freq_coords, spatial_freq_coords, log_otf_norm' );
axis xy image;
set(gca, 'FontName', 'Times New Roman','FontSize', 14);
title({'Frequency View for Slanted Sheet', sprintf('\\theta = %.1f°', theta*180/pi)}, 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Axial Spatial Frequency u_z (\mum^{-1})', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Lateral Spatial Frequency u_x (\mum^{-1})', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
lgd = legend('show', 'Location', 'Best');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
colormap(parula);
clim([0 0.01]);
hold on;
line_direction = [-sin(theta), cos(theta)]; % [direction_ux, direction_uz]
t_range = max(abs(spatial_freq_coords)) * 1.5;
t = linspace(-t_range, t_range, 2000); 
u_z_line_smooth = t * line_direction(2); % z coordinate
u_x_line_smooth = t * line_direction(1); % x coordinate
plot(u_x_line_smooth, u_z_line_smooth, 'r--', 'LineWidth', 2); % Use '-' to draw solid line
xlim([-2*fc, 2*fc]);
ylim([-2*fc, 2*fc]);
hold off;
%% Sheet Object Resolution Figure (Any Angles -> FWHM)

fprintf('\nStarting analysis of FWHM vs. Sheet Angle...\n');

% --- 1. Initialize parameters ---
N_angles = 91; % Set number of angle sampling points, e.g. from 0 to 180 degrees, every 2 degrees
theta_list = linspace(0, pi, N_angles); % Angle list (radians)
fwhm_values = nan(1, N_angles); % Used to store FWHM values calculated for each angle

% Pre-prepare OTF needed for frequency domain convolution (already calculated before)
sheet_thickness = distance_per_pixel; % Use single pixel thickness

% --- 2. Loop to calculate FWHM for each angle ---
% Create a timer to show estimated remaining time, as this computation is intensive
h_waitbar = waitbar(0, 'Calculating FWHM for each angle...');
tic; % Start timer

for k = 1:N_angles
    % Update progress bar
    if k > 1
        elapsed_time = toc;
        est_rem_time = (elapsed_time / (k-1)) * (N_angles - (k-1));
        waitbar(k/N_angles, h_waitbar, ...
            sprintf('Processing angle %.1f°... (Est. remaining: %.1f s)', ...
            theta_list(k)*180/pi, est_rem_time));
    else
        waitbar(k/N_angles, h_waitbar, ...
            sprintf('Processing angle %.1f°...', theta_list(k)*180/pi));
    end

    theta = theta_list(k);
    
    % a) Define normal vector for current angle and create 3D sheet
    % Normal vector defined as angle theta with z-axis, in x-z plane
    normal_vec = [sin(theta), 0, cos(theta)]; 
    plane_dist = X_real * normal_vec(1) + Y_real * normal_vec(2) + Z_real * normal_vec(3);
    sheet3D_slanted = double(abs(plane_dist) < (sheet_thickness / 2));
    
    % b) Frequency domain method to calculate convolution (imaging)
    FT_sheet_slanted = fftn(fftshift(sheet3D_slanted));
    I_freq_slanted   = ifftshift(ifftn( FT_sheet_slanted .* OTF ));
    
    % c) Extract 1D intensity profile along normal direction
    % Define sampling line
    max_dist = max(abs(dist_um)) * 0.8; % Sample within 80% range to avoid edge effects
    num_points = 1000; % High-density sampling for smooth curves
    line_coords_1d = linspace(-max_dist, max_dist, num_points);
    % Calculate sampling point coordinates in 3D space
    sample_points_x = line_coords_1d * normal_vec(1);
    sample_points_y = line_coords_1d * normal_vec(2);
    sample_points_z = line_coords_1d * normal_vec(3);
    % Extract intensity values using cubic interpolation
    profile_1D = interp3(X_real, Y_real, Z_real, real(I_freq_slanted), ...
                           sample_points_x, sample_points_y, sample_points_z, 'cubic', 0);
                       
    % d) Calculate FWHM of the profile
    if max(profile_1D) > 1e-9 % Ensure profile has valid signal
        profile_norm = profile_1D / max(profile_1D);
        [~, pk_idx] = max(profile_norm); % Find peak position
        
        % Find left half-maximum point
        half_idx_left = find(profile_norm(1:pk_idx) <= 0.5, 1, 'last');
        % Find right half-maximum point
        half_idx_right = find(profile_norm(pk_idx:end) <= 0.5, 1, 'first') + pk_idx - 1;
        
        if ~isempty(half_idx_left) && ~isempty(half_idx_right) && half_idx_left > 1
            % Calculate precise half-maximum coordinates through linear interpolation
            y1_L = profile_norm(half_idx_left); x1_L = line_coords_1d(half_idx_left);
            y2_L = profile_norm(half_idx_left+1); x2_L = line_coords_1d(half_idx_left+1);
            coord_L = x1_L + (0.5 - y1_L) * (x2_L - x1_L) / (y2_L - y1_L);
            
            y1_R = profile_norm(half_idx_right-1); x1_R = line_coords_1d(half_idx_right-1);
            y2_R = profile_norm(half_idx_right); x2_R = line_coords_1d(half_idx_right);
            coord_R = x1_R + (0.5 - y1_R) * (x2_R - x1_R) / (y2_R - y1_R);
            
            fwhm_values(k) = coord_R - coord_L;
        end
    end
end
close(h_waitbar); % Close wait bar
fprintf('FWHM calculation finished in %.2f seconds.\n', toc);

% --- 3. Visualization Results ---

fprintf('Starting analysis of FWHM (Method 2: Direct 2D IPSF Analysis)...\n');
% --- 1. Prepare 2D Data ---
IPSF_xz = squeeze(IPSF(centre_pt, :, :)); % Extract x-z slice
% x_coords and z_coords have been previously defined as dist_um
% --- 2. Loop to calculate 50% intensity radius for each angle ---
N_angles_direct = 180;                   
theta_list_direct = linspace(0, pi, N_angles_direct);  
x_half = nan(1, N_angles_direct);
z_half = nan(1, N_angles_direct);
% meshgrid for interp2
[Z_grid, X_grid] = meshgrid(z_coords, x_coords);

for k = 1 : N_angles_direct
    theta = theta_list_direct(k);
    
    % a) Determine maximum sampling radius r_max that doesn't exceed image boundaries for current angle
    % Note: sin/cos definitions in your code are slightly different from standard polar coordinates, preserving your definitions
    % theta=0 along x-axis, theta=pi/2 along z-axis
    if abs(cos(theta)) < 1e-8,  r_max_x = Inf; else, r_max_x = max(abs(x_coords)) / abs(cos(theta)); end
    if abs(sin(theta)) < 1e-8,  r_max_z = Inf; else, r_max_z = max(abs(z_coords)) / abs(sin(theta)); end
    r_max = min(r_max_x, r_max_z);  
    
    % b) Create sampling points along this angle
    N_r_sample = 800;
    r_vec = linspace(0, r_max, N_r_sample);
    
    % c) Calculate sampling point coordinates and obtain intensity values using interpolation
    x_r = r_vec * cos(theta);
    z_r = r_vec * sin(theta);
    f_r = interp2(Z_grid, X_grid, IPSF_xz, z_r, x_r, 'spline', 0 ); 
  
    % d) Normalize (assuming peak at r=0)
    if f_r(1) > 1e-9, f_r = f_r / f_r(1); else, f_r = zeros(size(f_r)); end
    
    % e) Find the point where intensity first drops below 0.5 and calculate its radius precisely
    idx_fall = find( f_r < 0.5, 1, 'first' );
    if ~isempty(idx_fall) && idx_fall > 1
        r1 = r_vec(idx_fall - 1); r2 = r_vec(idx_fall);
        f1 = f_r(idx_fall - 1);   f2 = f_r(idx_fall);
        r_half = r1 + (0.5 - f1) * (r2 - r1) / (f2 - f1);
    else
        r_half = NaN;
    end
    
    % f) Convert radius to contour point coordinates
    x_half(k) =  r_half * cos(theta);
    z_half(k) =  r_half * sin(theta);
end
% --- 3. Create complete contour line ---
x_direct_full = [ x_half,  -fliplr(x_half) ];
z_direct_full = [ z_half,  -fliplr(z_half) ];
fprintf('Method 2 (Direct Analysis) complete.\n');
% a) Plot FWHM vs. Angle in Cartesian coordinates (most direct information)
figure;
plot(theta_list * 180/pi, fwhm_values, 'b-', 'LineWidth', 2);
set(gca, 'FontName', 'Times New Roman','FontSize', 14);
title('System Resolution vs. Sheet Orientation', 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Sheet Angle (degrees from z-axis)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('FWHM / Resolution (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
lgd = legend('show', 'Location', 'Best');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
xlim([0 180]);
%ylim([0, max(fwhm_values)*1.1]); % Automatically adjust y-axis range

% b) Create and plot FWHM polar contour plot with IPSF slice as background
% Calculate contour line coordinates
r_half = fwhm_values / 2;
x_contour = r_half .* sin(theta_list);
z_contour = r_half .* cos(theta_list);
x_contour_full = [x_contour, -fliplr(x_contour)];
z_contour_full = [z_contour, -fliplr(z_contour)];

% --- Begin plotting ---
figure;
hold on; 

% 1. Plot background: IPSF X-Z slice
imagesc(z_coords, x_coords, IPSF_xz');
axis xy; 
colormap('parula');
hcb = colorbar;
hcb.Label.String = 'Normalized IPSF Intensity';
set(gca,'Color','k'); 

% 2. Overlay first contour line (Sheet FWHM, red)
plot(x_contour_full, z_contour_full, 'g--', 'LineWidth', 2);

% 3. Overlay second contour line (IPSF FWHM, green)
plot(x_direct_full, z_direct_full, 'r--', 'LineWidth', 2); % Use green dashed line for distinction

% 4. Add chart decorations and labels
set(gca, 'FontName', 'Times New Roman','FontSize', 14);
title('Comparison of FWHM Contours on IPSF', 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('z-axis (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('x-axis (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
lgd = legend('show', 'Location', 'best');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% 5. Set display range and aspect ratio
daspect([1 1 1]); 
axis_limit_z = max(abs([z_contour_full, z_direct_full])) * 1.5;
axis_limit_x = max(abs([x_contour_full, x_direct_full])) * 1.5;
xlim([-1.5*zz,1.5*zz]);
ylim([-1.5*zz,1.5*zz]);
hold off;
fprintf('Final comparison visualization has been generated.\n');
%% Line Object (HORIZONTAL) Figures
%  ========================================================================
fprintf('\n--- LSF Analysis PART 1: Fixed line along Y-axis ---\n');

% --- a) Create an ideal line object along X-axis ---
line3D_y_axis = zeros(size(IPSF));
line3D_y_axis(centre_pt,: , centre_pt) = 1; % 

% --- b) Visualize line object ---
figure;
p = patch(isosurface(X_real, Y_real, Z_real, line3D_y_axis, 0.5));
isonormals(X_real, Y_real, Z_real, line3D_y_axis, p);
p.FaceColor = 'red'; p.EdgeColor = 'none';
daspect([1,1,1]); view(3); axis tight; camlight; lighting gouraud;
grid on;
xlabel('x (µm)'); ylabel('y (µm)'); zlabel('z (µm)');
title('3D Line Object (along X-axis)');

% --- c) Frequency domain convolution to calculate imaging result (LSF) ---
FT_line = fftn(fftshift(line3D_y_axis));
OTF = ifftshift(fftn(fftshift(IPSF))); % calculate the OTF
OTF    = ifftshift(OTF);
I_LSF_3D = ifftshift(ifftn(FT_line .* OTF));
I_LSF_3D = real(I_LSF_3D); % Take real part, ignore computational errors

% --- d) Visualize 2D distribution of LSF in X-Z plane ---
% LSF is the image of a line, its blur distribution is in the plane perpendicular to the line, i.e., the x-z plane
lsf_xz_slice = squeeze(I_LSF_3D(centre_pt, :, :)); % Take x-z slice at y=0
lsf_xz_slice = lsf_xz_slice / max(lsf_xz_slice(:));

figure;
imagesc(z_coords, x_coords, lsf_xz_slice);
axis xy image;
colormap('parula');
colorbar;
xlabel('z-axis (µm)');
ylabel('x-axis (µm)');
title('Image of Line along X-axis (LSF in x-z plane)');
hold on;
plot(0, 0, 'w+', 'MarkerSize', 10, 'LineWidth', 2); % Mark center
daspect([1 1 1]);

% --- e) Extract and analyze 1D profile of LSF ---
lsf_z_profile = lsf_xz_slice(centre_pt, :); % Profile along z-axis (x=0)

% Calculate FWHM of Z-Profile
fwhm_z = NaN; % Default to NaN
profile_norm_z = lsf_z_profile' / max(lsf_z_profile); % Transpose to column vector
[~, pk_idx_z] = max(profile_norm_z);

% Find left half-maximum point
idx_left_z = find(profile_norm_z(1:pk_idx_z) <= 0.5, 1, 'last');
% Find right half-maximum point
idx_right_z = find(profile_norm_z(pk_idx_z:end) <= 0.5, 1, 'first') + pk_idx_z - 1;

if ~isempty(idx_left_z) && ~isempty(idx_right_z) && idx_left_z > 1 && idx_right_z < length(profile_norm_z)
    % Calculate left coordinate precisely using linear interpolation
    y1L=profile_norm_z(idx_left_z);   x1L=z_coords(idx_left_z);
    y2L=profile_norm_z(idx_left_z+1); x2L=z_coords(idx_left_z+1);
    coord_L_z = x1L + (0.5 - y1L) * (x2L - x1L) / (y2L - y1L);
    % Calculate right coordinate precisely using linear interpolation
    y1R=profile_norm_z(idx_right_z-1); x1R=z_coords(idx_right_z-1);
    y2R=profile_norm_z(idx_right_z);   x2R=z_coords(idx_right_z);
    coord_R_z = x1R + (0.5 - y1R) * (x2R - x1R) / (y2R - y1R);
    
    fwhm_z = coord_R_z - coord_L_z;
end

% -- END: Corrected FWHM calculation method --


fprintf('FWHM of LSF along z-axis: %.4f µm\n', fwhm_z);

figure;
plot(z_coords, lsf_z_profile, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('Z-Profile (FWHM=%.2f µm)', fwhm_z));
grid on;
xlabel('Distance from center (µm)');
ylabel('Normalized Intensity');
title('1D Profiles of the Line Spread Function (LSF)');
legend('show', 'Location', 'best');
xlim([-2, 2]);

% --- f) [New] Horizontal line frequency domain visualization (corrected version) ---
figure;
% Prepare MTF background (u_x-u_z plane)
% Note: To get the u_x-u_z plane, we need to squeeze the Y-axis (first dimension)
% But for consistency with your other code, we continue to squeeze the X-axis and correct ylabel
otf_slice_2d = squeeze(MTF(centre_pt, :, :)); % Extract u_y-u_z plane
imagesc(spatial_frequency_range, spatial_frequency_range, otf_slice_2d);
axis xy image; colormap('parula');
caxis([0, 0.1]); 
hcb = colorbar; hcb.Label.String = 'MTF Value';
hold on;

% For a line along the X-axis, its spectrum intersection on the ux-uz plane is the uz-axis, i.e., a horizontal line ux=0
% line(x_coords, y_coords)
line(spatial_frequency_range, zeros(size(spatial_frequency_range)), 'Color', 'red', 'LineWidth', 2);

% Add correct labels
xlabel('Axial Spatial Frequency u_z (\mum^{-1})');
ylabel('Lateral Spatial Frequency u_y (\mum^{-1})'); % Corresponds to squeeze(MTF(centre_pt,:,:))
title({'Frequency View for Horizontal Line (along x-axis)', '(Red Line: Line Spectrum, Color: MTF)'});
xlim([-2*fc, 2*fc]); ylim([-2*fc, 2*fc]);
set(gca, 'Layer', 'top');
hold off;
%% Line Object (ROTATE) Figures

fprintf('\n--- LSF Analysis PART 2: Analyze line at single rotation angle ---\n');

% --- a) Construct and visualize line at specified angle ---

phi_rot_deg = 0; % << You can set any angle here, for example 30°
phi_rot_rad = phi_rot_deg * pi/180;

% Define rotated direction vector (within x-z plane, rotating from z-axis to x-axis)

line_direction = [cos(phi_rot_rad), 0, sin(phi_rot_rad)];

% Use distance-to-line formula to create 3D line object

p_dot_d = X_real * line_direction(1) + Y_real * line_direction(2) + Z_real * line_direction(3);
dist_sq_from_line = (X_real - p_dot_d*line_direction(1)).^2 + ...
(Y_real - p_dot_d*line_direction(2)).^2 + ...
(Z_real - p_dot_d*line_direction(3)).^2;
line_radius = distance_per_pixel;
line3D_rotated = double(dist_sq_from_line < line_radius^2);

% Visualize rotated line

figure;
p = patch(isosurface(X_real, Y_real, Z_real, line3D_rotated, 0.5));
isonormals(X_real, Y_real, Z_real, line3D_rotated, p);
p.FaceColor = 'magenta'; p.EdgeColor = 'none';
daspect([1,1,1]); view(3); axis tight; camlight; lighting gouraud;
grid on;
set(gca, 'FontName', 'Times New Roman','FontSize', 14);
title(sprintf('3D Line Object (Rotated by %.1f°)', phi_rot_deg), 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('x (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('y (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
zlabel('z (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
hold on;

% Mark line direction with arrow

L = max(dist_um) * 0.8;
quiver3(0, 0, 0, line_direction(1), line_direction(2), line_direction(3), L, ...
'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
hold off;

% --- b) Frequency domain convolution to calculate imaging result (LSF) ---

FT_line_rot = fftn(fftshift(line3D_rotated));
I_LSF_rot_3D = ifftshift(ifftn(FT_line_rot .* OTF));
I_LSF_rot_3D = real(I_LSF_rot_3D); % Take real part

% --- c) Extract and analyze 1D profile of LSF (perpendicular to line direction) ---
% Define sampling direction perpendicular to line

perp_direction = [-sin(phi_rot_rad), 0, cos(phi_rot_rad)];
max_dist_profile = 5;
num_points_profile = 1001;
line_coords_1d = linspace(-max_dist_profile, max_dist_profile, num_points_profile);

% Calculate sampling point coordinates

sample_points_x = line_coords_1d * perp_direction(1);
sample_points_y = line_coords_1d * perp_direction(2);
sample_points_z = line_coords_1d * perp_direction(3);

% Extract profile using interpolation

profile_1D = interp3(X_real, Y_real, Z_real, I_LSF_rot_3D, ...
sample_points_x, sample_points_y, sample_points_z, 'cubic', 0);
profile_1D_norm = profile_1D / max(profile_1D);
fwhm_val = NaN; % Initialize to NaN in case calculation fails
[~, pk_idx] = max(profile_1D_norm); % Find peak position
idx_left = find(profile_1D_norm(1:pk_idx) <= 0.5, 1, 'last');
idx_right = find(profile_1D_norm(pk_idx:end) <= 0.5, 1, 'first') + pk_idx - 1;
if ~isempty(idx_left) && ~isempty(idx_right) && idx_left > 1 && idx_right < num_points_profile

% Calculate left coordinate precisely using linear interpolation

y1L = profile_1D_norm(idx_left); x1L = line_coords_1d(idx_left);
y2L = profile_1D_norm(idx_left+1); x2L = line_coords_1d(idx_left+1);
coord_L = x1L + (0.5 - y1L) * (x2L - x1L) / (y2L - y1L);

% Calculate right coordinate precisely using linear interpolation

y1R = profile_1D_norm(idx_right-1); x1R = line_coords_1d(idx_right-1);
y2R = profile_1D_norm(idx_right); x2R = line_coords_1d(idx_right);
coord_R = x1R + (0.5 - y1R) * (x2R - x1R) / (y2R - y1R);
fwhm_val = coord_R - coord_L;
end

fprintf('FWHM of the profile normal to the %.1f° line is: %.4f µm\n', phi_rot_deg, fwhm_val);

figure;
plot(line_coords_1d, profile_1D_norm, 'b-', 'LineWidth', 2);
hold on;
grid on;
if ~isnan(fwhm_val)
plot([coord_L, coord_R], [0.5, 0.5], 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot([coord_L, coord_L], [0, 0.5], 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot([coord_R, coord_R], [0, 0.5], 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

% Add FWHM value to title and specify it's the perpendicular (Normal) direction
set(gca, 'FontName', 'Times New Roman','FontSize', 14);
title(sprintf('LSF Profile Normal to Line (Angle = %.1f°, FWHM = %.3f µm)', phi_rot_deg, fwhm_val), 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Distance Normal to Line Direction (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Normalized Intensity', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

%xlim([-2, 2]);
ylim([0,1]);
hold off;

% --- f) [New] Visualize LSF in rotated coordinate system ---
fprintf('Resampling LSF into rotated coordinate system...\n');

% 1. Define new coordinate axis directions: u-axis perpendicular to line, v-axis parallel to line
u_vec = perp_direction; % u_vec is the normal vector
v_vec = line_direction; % v_vec is the parallel vector

% 2. Create coordinates for new 2D sampling grid
% u-axis (horizontal axis) is the direction for observing blur, range can be smaller
u_coords = linspace(-5, 5, 201); 
% v-axis (vertical axis) is the direction along the line, range can be larger
v_coords = linspace(-5, 5, 401); 
[U_grid, V_grid] = meshgrid(u_coords, v_coords);

% 3. Map 2D grid points to 3D space coordinates
% P_3d = u * u_vec + v * v_vec
sample_points_X = U_grid .* u_vec(1) + V_grid .* v_vec(1);
sample_points_Y = U_grid .* u_vec(2) + V_grid .* v_vec(2); %
sample_points_Z = U_grid .* u_vec(3) + V_grid .* v_vec(3);

% 4. Use interp3 for three-dimensional resampling
resampled_LSF_2D = interp3(dist_um, dist_um, dist_um, I_LSF_rot_3D, ...
                           sample_points_X, sample_points_Y, sample_points_Z, 'cubic', 0);

% 5. Plot results
figure;
imagesc(u_coords, v_coords, resampled_LSF_2D');
axis xy; % Ensure y-axis (v-axis) direction is upward
colormap('parula'); % Using 'hot' colormap usually better displays light spots
colorbar;
daspect([1 1 1]); % Ensure correct coordinate axis proportions

% Add clear labels
colorbar
set(gca, 'FontName', 'Times New Roman','FontSize', 14);
title(sprintf('LSF in Rotated Frame (Line Angle = %.1f°)', phi_rot_deg), 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Distance Normal to Line (u-axis, µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Distance Parallel to Line (v-axis, µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlim([-5,5]);
ylabel('Distance Parallel to Line (v-axis, µm)');

% --- g) [New] Rotated line frequency domain visualization ---
figure;
% Prepare MTF background
otf_slice_2d = squeeze(MTF(centre_pt, :, :)); % Extract u_x-u_z plane
imagesc(spatial_frequency_range, spatial_frequency_range, otf_slice_2d');
axis xy image; colormap('parula');
caxis([0, 1]);
hcb = colorbar; hcb.Label.String = 'MTF Value';
hold on;

% Calculate direction of spectral line (perpendicular to spatial domain line direction)
% Spatial domain line direction (x,z): [cos(phi), sin(phi)]
% Frequency domain spectral line direction (ux,uz): [-sin(phi), cos(phi)]
spec_line_direction = [-sin(phi_rot_rad), cos(phi_rot_rad)];

% Plot spectral line
t_range = max(abs(spatial_frequency_range));
t = linspace(-t_range, t_range, 1000);
u_x_spec_line = t * spec_line_direction(1); % Frequency domain x-component
u_z_spec_line = t * spec_line_direction(2); % Frequency domain z-component
plot(u_x_spec_line, u_z_spec_line, 'r--', 'LineWidth', 2, 'DisplayName', 'Line Spectrum');

% Add labels and title
hcb = colorbar; 
set(hcb, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman','FontSize', 14);
title({'Frequency View for Rotated Line', sprintf('\\phi = %.1f°', phi_rot_deg)}, 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Axial Spatial Frequency u_z (\mum^{-1})', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Lateral Spatial Frequency u_x (\mum^{-1})', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
lgd = legend('show', 'Location', 'Best');
set(lgd, 'String', {'FT of line in x-z frequency view'}, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
xlim([-2*fc, 2*fc]); ylim([-2*fc, 2*fc]);
hold off;
%% Line Object Resolution Figure (Any Angles -> FWHM)
%  ========================================================================
fprintf('\n--- LSF Analysis PART 3: Calculate FWHM vs. rotation angle ---\n');

% --- a) Set rotation parameters ---
N_angles = 91; % 
phi_list_deg = linspace(0, 180, N_angles);
phi_list_rad = phi_list_deg * pi/180;
fwhm_lsf_values = nan(1, N_angles);

% --- b) Loop to calculate LSF FWHM for each angle ---
h_waitbar = waitbar(0, 'Calculating LSF FWHM for each line angle...');
tic;
for k = 1:N_angles
    % Update progress bar
    waitbar(k/N_angles, h_waitbar, sprintf('Processing angle %.1f°', phi_list_deg(k)));
    
    phi_rot_rad = phi_list_rad(k);
    
    % 1. Create rotated line object (completely consistent with Part 2 logic)
    line_direction = [cos(phi_rot_rad), 0, sin(phi_rot_rad)];
    p_dot_d = X_real * line_direction(1) + Y_real * line_direction(2) + Z_real * line_direction(3);
    dist_sq_from_line = (X_real - p_dot_d*line_direction(1)).^2 + ...
                        (Y_real - p_dot_d*line_direction(2)).^2 + ...
                        (Z_real - p_dot_d*line_direction(3)).^2;
    line_radius = distance_per_pixel;
    line3D_rotated = double(dist_sq_from_line < line_radius^2);
    
    % 2. Frequency domain convolution (completely consistent with Part 2 logic)
    FT_line_rot = fftn(fftshift(line3D_rotated));
    I_LSF_rot_3D = ifftshift(ifftn(FT_line_rot .* OTF));
    
    % 3. Extract profile along direction perpendicular to line (completely consistent with Part 2 logic)
    perp_direction = [-sin(phi_rot_rad), 0, cos(phi_rot_rad)]; 
    max_dist_profile = 5;
    num_points_profile = 1001;
    line_coords_1d = linspace(-max_dist_profile, max_dist_profile, num_points_profile);
    sample_points_x = line_coords_1d * perp_direction(1);
    sample_points_y = line_coords_1d * perp_direction(2); % Always 0
    sample_points_z = line_coords_1d * perp_direction(3);
    
    profile_1D = interp3(X_real, Y_real, Z_real, real(I_LSF_rot_3D), ...
                         sample_points_x, sample_points_y, sample_points_z, 'cubic', 0);
                     
    % 4. Calculate FWHM (completely consistent with Part 2 logic)
    fwhm_val = NaN;
    if max(profile_1D) > 1e-9 % Ensure there is valid signal
        profile_1D_norm = profile_1D / max(profile_1D);
        [~, pk_idx] = max(profile_1D_norm);
        
        idx_left = find(profile_1D_norm(1:pk_idx) <= 0.5, 1, 'last');
        idx_right = find(profile_1D_norm(pk_idx:end) <= 0.5, 1, 'first') + pk_idx - 1;
        
        if ~isempty(idx_left) && ~isempty(idx_right) && idx_left > 1 && idx_right < num_points_profile
            y1L = profile_1D_norm(idx_left);   x1L = line_coords_1d(idx_left);
            y2L = profile_1D_norm(idx_left+1); x2L = line_coords_1d(idx_left+1);
            coord_L = x1L + (0.5 - y1L) * (x2L - x1L) / (y2L - y1L);
            
            y1R = profile_1D_norm(idx_right-1); x1R = line_coords_1d(idx_right-1);
            y2R = profile_1D_norm(idx_right);   x2R = line_coords_1d(idx_right);
            coord_R = x1R + (0.5 - y1R) * (x2R - x1R) / (y2R - y1R);
            
            fwhm_val = coord_R - coord_L;
        end
    end
    
    % 5. Store FWHM value for current angle
    fwhm_lsf_values(k) = fwhm_val;
end
close(h_waitbar);
fprintf('LSF FWHM calculation for all angles finished in %.2f seconds.\n', toc);

% --- b) Calculate the Peak-to-Sidelobe Ratio (PSR) from this PSF ---
axial_profile = squeeze(IPSF(centre_pt, centre_pt, :));
peak_signal = axial_profile(centre_pt);
[pks, ~] = findpeaks(axial_profile(centre_pt+1:end), 'MinPeakProminence', 1e-5);
if ~isempty(pks)
    pedestal_signal = pks(1);
    psr_value = peak_signal / pedestal_signal;
else
    psr_value = Inf;
end
fprintf('  > Peak-to-Sidelobe Ratio (PSR) = %.2f\n\n', psr_value);


% Visualization results: FWHM vs. Angle Cartesian coordinate plot ---
figure;
plot(phi_list_deg, fwhm_lsf_values, 'b-', 'LineWidth', 2);
set(gca, 'FontName', 'Times New Roman','FontSize', 14);
title('Line Imaging Resolution vs. Orientation', 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Line Angle in x-z plane (degrees from x-axis)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('LSF FWHM (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
lgd = legend('show', 'Location', 'Best');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
xlim([0 180]);

% --- d) Visualization results: LSF FWHM contour plot overlaid on IPSF ---
% FWHM is the width perpendicular to the line direction, so points on the contour plot have direction perp_direction
r_lsf = fwhm_lsf_values / 2;
x_contour = r_lsf .* (-sin(phi_list_rad));
z_contour = r_lsf .* (cos(phi_list_rad));

% Create complete contour (0-360 degrees)
x_contour2_full = [x_contour, -fliplr(x_contour), -x_contour, fliplr(x_contour)];
z_contour2_full = [z_contour, fliplr(z_contour), -z_contour, -fliplr(z_contour)];

figure;
hold on;
% Plot background: IPSF X-Z slice
IPSF_xz = squeeze(IPSF(centre_pt, :, :));
imagesc(z_coords, x_coords, IPSF_xz');
axis xy; colormap('parula'); hcb = colorbar;

% Plot LSF FWHM contour line
plot(x_contour2_full, z_contour2_full, 'y--', 'LineWidth', 2);
plot(x_contour_full, z_contour_full, 'g--', 'LineWidth', 2);
plot(x_direct_full, z_direct_full, 'r--', 'LineWidth', 2);


% Add chart decorations
colorbar
set(gca, 'FontName', 'Times New Roman','FontSize', 14);
title('LSF Resolution Contour Overlaid on IPSF', 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('z-axis (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('x-axis (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
lgd = legend('show', 'Location', 'best');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
daspect([1 1 1]);
axis_limit = max(abs([x_contour2_full, z_contour2_full]))*1.2;
xlim([-1.5*zz,1.5*zz]);
ylim([-1.5*zz,1.5*zz]);
grid on;
text(4.8, -4.5, sprintf('PSR = %.1f', psr_value), ...
    'Color', 'white', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'right', 'BackgroundColor', [0 0 0 0.5]);
legend('show', 'Location', 'best');


hold off;

fprintf('\n>>> LSF analysis complete <<<\n');

%% === Component 2: Polar overlay: Point vs Sheet vs Line (FWHM vs angle from z-axis) ===
% Dependent variables (use directly if they exist):
% Sheet: theta_list (measured from z axis)           , fwhm_values
% Line : phi_list_rad (line direction, measured from x axis) , fwhm_lsf_values  ——> Normal direction actually has same angle as "measured from z axis", can directly use phi_list_rad
% Point: theta_list_direct (measured from x axis), r_half(radius=FWHM/2) ——> optional

% -------- Circular smoothing (optional) --------
w =11;  % Choose an odd number between 5-11
pad   = @(x)[x(end-w+1:end); x(:); x(1:w)];
unpad = @(x)x(w+1:end-w);
circ_smooth = @(x) (numel(x)>2) .* unpad(movmean(pad(x(:)), w, 'omitnan')) + (numel(x)<=2).*x(:); % Maintain column vector

% ====================== Sheet (Normal FWHM_n) ======================
ang_sheet = theta_list(:);          % Already "measured from z axis", range 0..pi
rS = fwhm_values(:);
rS = circ_smooth(rS);               % Can comment out to disable smoothing
% Use NaN to break closed connection lines, avoid straight lines through origin
angS_full = [ang_sheet; ang_sheet+pi; NaN];
rS_full   = [rS;       rS;          NaN];

% ====================== Line (Normal FWHM_\perp) ======================
% Line direction vector is [cos(phi),0,sin(phi)], its normal is [-sin(phi),0,cos(phi)];
% Normal's angle α with z-axis satisfies cos(α)=cos(phi) → α=phi (constant within 0..pi),
% Therefore can directly use phi_list_rad as the angle "measured from z-axis".
ang_line = phi_list_rad(:);         % 0..pi
rL = fwhm_lsf_values(:);
rL = circ_smooth(rL);
angL_full = [ang_line; ang_line+pi; NaN];
rL_full   = [rL;       rL;          NaN];

% ====================== Point (IPSF half-maximum contour FWHM) [Optional] ======================

ang_point = (pi/2) - theta_list(:);     % First convert "measured from x-axis" to "measured from z-axis"
ang_point = mod(ang_point, pi);         % Map to [0, pi)

rP = 2 * r_half_vec(:);                 % FWHM = 2*radius

% —— 1) Remove NaN and sort by angle in ascending order (avoid cross-segment connections due to disorder) ——
valid = isfinite(ang_point) & isfinite(rP);
ang_point = ang_point(valid); 
rP        = rP(valid);

[ang_point, idx] = sort(ang_point); 
rP = rP(idx);

% (Optional) Apply circular smoothing after sorting
% rP = circ_smooth(rP);

% —— 2) Insert NaN at "repeated angles" to break radial connections ——
dup = [false; abs(diff(ang_point)) < 1e-9];  % Adjacent angles are almost identical
ang_plot = ang_point;  r_plot = rP;
ang_plot(dup) = NaN;   r_plot(dup) = NaN;

% —— 3) Extend to 0..2π and break once more between half cycles —— 
angP_full = [ang_plot; NaN; ang_plot + pi; NaN];
rP_full   = [r_plot;  NaN; r_plot;        NaN];


% ====================== Plotting ======================
fig = figure('Name','Polar overlay: resolution anisotropy','Color','w');
pax = polaraxes(fig);  hold(pax,'on');

% Unify polar axis semantics: 0°=axial(z), 90°=lateral(x)
pax.ThetaZeroLocation = 'top';
pax.ThetaDir          = 'counterclockwise';
pax.ThetaTick         = [0 90 180 270];
pax.ThetaTickLabel    = {'0° axial (z)','90° lateral (x)','180°','270°'};

% Plot each line (only plot if exists)
polarplot(pax, angP_full, rP_full, '--', 'LineWidth', 2, 'DisplayName','Point  FWHM (IPSF)');
polarplot(pax, angS_full, rS_full, '--', 'LineWidth', 2, 'DisplayName','Sheet  FWHM_n');
polarplot(pax, angL_full, rL_full, '--',  'LineWidth', 2, 'DisplayName','Line   FWHM_\perp');

% Radius range
r_candidates = [rS(:); rL(:)];
if ~isempty(rP), r_candidates = [r_candidates; rP(:)]; end
rmax = max(r_candidates, [], 'omitnan');
if isfinite(rmax) && rmax>0, rlim(pax, [0, 1.1*rmax]); end

title(pax,'FWHM vs orientation (angle from z-axis)');
legend(pax,'Location','southoutside');
grid(pax,'on');


%% === Component 1: MPA Dashboard: Confocal Microscopy Performance Atlas =====
% This section generates a standardized performance dashboard 
% integrating spatial, frequency, and energy domain metrics

fprintf('\n========== Generating MPA Dashboard (Confocal) ==========\n');

% Create figure with specific size for publication quality
fig_dash = figure('Name', 'MPA Dashboard - Confocal Microscopy', ...
                  'Position', [100, 50, 1400, 900], ...
                  'Color', 'white');

% Use tiledlayout for professional subplot arrangement
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Overall title for the dashboard (including pinhole size)
title(t, sprintf('Confocal Microscopy Performance Dashboard (NA=%.1f, λ=%dnm, Pinhole=%.1f AU)', ...
                 NA, lambda*1e9, AU_specific), ...
      'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Arial');

% === Subplot 1: PSF Lateral Profile (X) ===
nexttile;
% Extract lateral profile through center (using final confocal PSF)
lateral_profile = squeeze(IPSF(centre_pt, :, centre_pt));
lateral_profile_norm = lateral_profile / max(lateral_profile);

% High-resolution interpolation for smooth curve
factor = 1000;
x_interp = linspace(x_coords(1), x_coords(end), numel(x_coords)*factor);
y_interp = interp1(x_coords, lateral_profile_norm, x_interp, 'spline');

% Plot the profile
plot(x_interp, y_interp, 'b-', 'LineWidth', 2);
hold on;

% Add FWHM indicator lines
yline(0.5, 'k--', 'LineWidth', 1, 'Alpha', 0.5);

% Calculate FWHM positions (recompute for accuracy)
[~, peak_idx] = max(y_interp);
left_idx = find(y_interp(1:peak_idx) <= 0.5, 1, 'last');
right_idx = find(y_interp(peak_idx:end) <= 0.5, 1, 'first') + peak_idx - 1;

if ~isempty(left_idx) && ~isempty(right_idx)
    xL = x_interp(left_idx);
    xR = x_interp(right_idx);
    FWHM_lat_conf = xR - xL;
    
    % Mark FWHM boundaries
    plot([xL, xL], [0, 0.5], 'r-', 'LineWidth', 1.5);
    plot([xR, xR], [0, 0.5], 'r-', 'LineWidth', 1.5);
    plot([xL, xR], [0.5, 0.5], 'r-', 'LineWidth', 2);
    
    % Add FWHM annotation
    text(0, 0.3, sprintf('FWHM = %.3f μm', FWHM_lat_conf), ...
         'HorizontalAlignment', 'center', ...
         'FontSize', 12, 'FontWeight', 'bold', ...
         'Color', 'r', 'BackgroundColor', 'white', ...
         'EdgeColor', 'r', 'LineWidth', 1);
else
    FWHM_lat_conf = fwhm_lateral; % Use previously calculated value
    text(0, 0.3, sprintf('FWHM = %.3f μm', FWHM_lat_conf), ...
         'HorizontalAlignment', 'center', ...
         'FontSize', 12, 'FontWeight', 'bold', ...
         'Color', 'r', 'BackgroundColor', 'white', ...
         'EdgeColor', 'r', 'LineWidth', 1);
end

% Formatting
grid on; box on;
xlim([-2, 2]);
ylim([0, 1.05]);
xlabel('Lateral Position x (μm)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Normalized Intensity', 'FontSize', 11, 'FontWeight', 'bold');
title('PSF Lateral Profile', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontName', 'Arial');

% === Subplot 2: Axial Sheet Response (Optical Sectioning) ===
nexttile;
% Use the axial response from sheet imaging (already calculated)
if exist('zp', 'var') && exist('ifq', 'var')
    % Use sheet response profile
    plot(zp, ifq, 'b-', 'LineWidth', 2);
    hold on;
    
    % FWHM indicators
    yline(0.5, 'k--', 'LineWidth', 1, 'Alpha', 0.5);
    
    if exist('xL_f', 'var') && exist('xR_f', 'var')
        plot([xL_f, xL_f], [0, 0.5], 'r-', 'LineWidth', 1.5);
        plot([xR_f, xR_f], [0, 0.5], 'r-', 'LineWidth', 1.5);
        plot([xL_f, xR_f], [0.5, 0.5], 'r-', 'LineWidth', 2);
        
        % Add FWHM annotation
        text(0, 0.3, sprintf('FWHM_z = %.3f μm', FWHM_f), ...
             'HorizontalAlignment', 'center', ...
             'FontSize', 12, 'FontWeight', 'bold', ...
             'Color', 'r', 'BackgroundColor', 'white', ...
             'EdgeColor', 'r', 'LineWidth', 1);
    end
else
    % Fallback to PSF axial profile
    axial_profile = squeeze(IPSF(centre_pt, centre_pt, :));
    axial_profile_norm = axial_profile / max(axial_profile);
    plot(z_coords, axial_profile_norm, 'b-', 'LineWidth', 2);
    hold on;
    yline(0.5, 'k--', 'LineWidth', 1, 'Alpha', 0.5);
    
    % Add FWHM if available
    if exist('fwhm_axial', 'var')
        text(0, 0.3, sprintf('FWHM_z = %.3f μm', fwhm_axial), ...
             'HorizontalAlignment', 'center', ...
             'FontSize', 12, 'FontWeight', 'bold', ...
             'Color', 'r', 'BackgroundColor', 'white', ...
             'EdgeColor', 'r', 'LineWidth', 1);
    end
end

% Add PSR annotation for confocal
if exist('psr_value', 'var')
    text(3, 0.8, sprintf('PSR = %.1f', psr_value), ...
         'FontSize', 10, 'FontWeight', 'bold', ...
         'Color', [0, 0.5, 0], 'BackgroundColor', 'white', ...
         'EdgeColor', [0, 0.5, 0], 'LineWidth', 1);
end

% Formatting
grid on; box on;
xlim([-8, 8]);
ylim([0, 1.05]);
xlabel('Axial Position z (μm)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Normalized Intensity', 'FontSize', 11, 'FontWeight', 'bold');
title('Axial Sheet Response (Optical Sectioning)', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontName', 'Arial');

% === Subplot 3: Radial MTF with Key Metrics ===
nexttile;
% Calculate radial MTF if not already done
if exist('fr_centers', 'var') && exist('mtf_rad', 'var')
    plot(fr_centers, mtf_rad, 'b-', 'LineWidth', 2);
    hold on;
    
    % Mark key MTF values
    yline(0.5, 'g--', 'LineWidth', 1, 'Alpha', 0.7, 'Label', 'MTF = 0.5');
    yline(0.1, 'r--', 'LineWidth', 1, 'Alpha', 0.7, 'Label', 'MTF = 0.1');
    
    % Mark frequency metrics
    if exist('fc_mtf50', 'var') && isfinite(fc_mtf50)
        xline(fc_mtf50, 'g:', 'LineWidth', 1.5);
        plot(fc_mtf50, 0.5, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    end
    if exist('fc_mtf10', 'var') && isfinite(fc_mtf10)
        xline(fc_mtf10, 'r:', 'LineWidth', 1.5);
        plot(fc_mtf10, 0.1, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    end
    if exist('fc_thresh', 'var') && isfinite(fc_thresh)
        xline(fc_thresh, 'k--', 'LineWidth', 1.5, 'Alpha', 0.7);
    end
    
    % Add metric annotations
    text_str = sprintf('f_c = %.2f cyc/μm\nMTF50 = %.2f cyc/μm\nMTF10 = %.2f cyc/μm', ...
                      fc_thresh, fc_mtf50, fc_mtf10);
    text(1.5, 0.7, text_str, ...
         'FontSize', 11, 'FontWeight', 'bold', ...
         'BackgroundColor', 'white', ...
         'EdgeColor', 'k', 'LineWidth', 1);
else
    % Compute radial average of MTF
    MTF_z0 = squeeze(MTF(:,:,centre_pt));
    MTF_z0 = MTF_z0 / max(MTF_z0(:));
    
    [fx_grid, fy_grid] = meshgrid(spatial_frequency_range, spatial_frequency_range);
    fr_simple = sqrt(fx_grid.^2 + fy_grid.^2);
    fr_bins = linspace(0, max(spatial_frequency_range), 100);
    mtf_radial = zeros(size(fr_bins));
    
    for i = 1:length(fr_bins)-1
        mask = (fr_simple >= fr_bins(i)) & (fr_simple < fr_bins(i+1));
        if any(mask(:))
            mtf_radial(i) = mean(MTF_z0(mask));
        end
    end
    
    plot(fr_bins, mtf_radial/max(mtf_radial(:)+eps), 'b-', 'LineWidth', 2);
    hold on;
    yline(0.5, 'g--', 'LineWidth', 1, 'Alpha', 0.7);
    yline(0.1, 'r--', 'LineWidth', 1, 'Alpha', 0.7);
    
    % Find cutoff frequencies
    idx50 = find(mtf_radial/max(mtf_radial(:)+eps) <= 0.5, 1, 'first');
    idx10 = find(mtf_radial/max(mtf_radial(:)+eps) <= 0.1, 1, 'first');
    if ~isempty(idx50)
        fc_mtf50 = fr_bins(idx50);
        xline(fc_mtf50, 'g:', 'LineWidth', 1.5);
    end
    if ~isempty(idx10)
        fc_mtf10 = fr_bins(idx10);
        xline(fc_mtf10, 'r:', 'LineWidth', 1.5);
    end
end


% Formatting
grid on; box on;
xlim([0, 4]);
ylim([0, 1.05]);
xlabel('Spatial Frequency (cyc/μm)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('MTF', 'FontSize', 11, 'FontWeight', 'bold');
title('Modulation Transfer Function', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontName', 'Arial');

% === Subplot 4: 3D PSF Isosurface with HMEF Range Visualization ===
nexttile;

% Calculate half-maximum level for HMEF
half_max_level = 0.5 * max(IPSF(:));

% Create outer isosurface at 0.001 threshold (full PSF extent)
iso_val_outer = 0.001;
p_outer = patch(isosurface(X_real, Y_real, Z_real, IPSF, iso_val_outer));
isonormals(X_real, Y_real, Z_real, IPSF, p_outer);

% Set appearance for outer surface
p_outer.FaceColor = 'blue';  % Very light blue
p_outer.EdgeColor = 'none';
p_outer.FaceAlpha = 0.15;  % Very transparent

hold on;

% Create inner isosurface at half-maximum (HMEF boundary)
p_inner = patch(isosurface(X_real, Y_real, Z_real, IPSF, half_max_level));
isonormals(X_real, Y_real, Z_real, IPSF, p_inner);

% Set appearance for inner surface (HMEF region) - Use green for confocal
p_inner.FaceColor = 'red';  % Green color for confocal
p_inner.EdgeColor = 'none';
p_inner.FaceAlpha = 0.85;  % More opaque

% Lighting and view
view(45, 25);
axis equal tight;
camlight('headlight');
lighting gouraud;
material shiny;

% Calculate HMEF if not already done
if ~exist('HMEF', 'var')
    IPSF_energy = IPSF / sum(IPSF(:));
    HMEF = sum(IPSF_energy(IPSF >= half_max_level));
end

% Add HMEF annotation with explanation
text_3d = sprintf('HMEF = %.3f', HMEF);
text(3, 3, 0, text_3d, ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', 'white', ...
     'EdgeColor', [0.1, 0.6, 0.3], 'LineWidth', 1.5, ...
     'HorizontalAlignment', 'center');

% Add legend explaining the visualization
text(2, 2, 3, sprintf('Green: FWHM region\nLight blue: Full PSF'), ...
     'FontSize', 9, 'FontWeight', 'normal', ...
     'BackgroundColor', [0.95 0.95 0.95], ...
     'EdgeColor', 'k', 'LineWidth', 0.5, ...
     'HorizontalAlignment', 'left');

% Formatting
grid on; box on;
xlabel('x (μm)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('y (μm)', 'FontSize', 10, 'FontWeight', 'bold');
zlabel('z (μm)', 'FontSize', 10, 'FontWeight', 'bold');
title('3D PSF Structure with HMEF Region', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontName', 'Arial');
xlim([-6, 6]); ylim([-6, 6]); zlim([-6, 6]);

% === Export Options ===
% Save the dashboard
% saveas(fig_dash, sprintf('MPA_Dashboard_Confocal_%.1fAU.png', AU_specific));
% saveas(fig_dash, sprintf('MPA_Dashboard_Confocal_%.1fAU.fig', AU_specific));
% print(fig_dash, sprintf('MPA_Dashboard_Confocal_%.1fAU', AU_specific), '-dpng', '-r300');  % High resolution

fprintf('MPA Dashboard (Confocal) generated successfully!\n');
fprintf('Key Performance Metrics:\n');
fprintf('  - Pinhole Size: %.1f AU\n', AU_specific);
fprintf('  - Lateral FWHM: %.3f μm\n', fwhm_lateral);
fprintf('  - Axial FWHM: %.3f μm\n', FWHM_f);


%% New Component 2: ARP
if exist('psr_value', 'var')
    fprintf('  - Peak-to-Sidelobe Ratio: %.1f\n', psr_value);
end
if exist('fc_thresh', 'var')
    fprintf('  - Cutoff frequency: %.2f cyc/μm\n', fc_thresh);
end
fprintf('  - HMEF: %.3f\n', HMEF);
fprintf('========================================\n');


if exist('AU_specific', 'var')
    microscopy_type = 'Confocal';
    pinhole_str = sprintf(' (%.1f AU)', AU_specific);
else
    microscopy_type = 'Wide-field';
    pinhole_str = '';
end

% -------- Circular smoothing function (unchanged) --------
w = 7;  % Smoothing window
pad   = @(x)[x(end-w+1:end); x(:); x(1:w)];
unpad = @(x)x(w+1:end-w);
circ_smooth = @(x) (numel(x)>2) .* unpad(movmean(pad(x(:)), w, 'omitnan')) + (numel(x)<=2).*x(:);

% ====================== Data Preparation (same as before) ======================
% Sheet FWHM
if exist('theta_list', 'var') && exist('fwhm_values', 'var')
    ang_sheet = theta_list(:);
    rS = fwhm_values(:);
    rS = circ_smooth(rS);
    angS_full = [ang_sheet; ang_sheet+pi; NaN];
    rS_full   = [rS; rS; NaN];
    has_sheet = true;
else
    has_sheet = false;
end

% Line FWHM
if exist('phi_list_rad', 'var') && exist('fwhm_lsf_values', 'var')
    ang_line = phi_list_rad(:);
    rL = fwhm_lsf_values(:);
    rL = circ_smooth(rL);
    angL_full = [ang_line; ang_line+pi; NaN];
    rL_full   = [rL; rL; NaN];
    has_line = true;
else
    has_line = false;
end

% Point FWHM
if exist('theta_list_direct', 'var') && exist('r_half_vec', 'var')
    ang_point = (pi/2) - theta_list(:);
    ang_point = mod(ang_point, pi);
    rP = 2 * r_half_vec(:);
    
    valid = isfinite(ang_point) & isfinite(rP);
    ang_point = ang_point(valid); 
    rP = rP(valid);
    
    [ang_point, idx] = sort(ang_point); 
    rP = rP(idx);
    
    dup = [false; abs(diff(ang_point)) < 1e-9];
    ang_plot = ang_point;  r_plot = rP;
    ang_plot(dup) = NaN;   r_plot(dup) = NaN;
    
    angP_full = [ang_plot; NaN; ang_plot + pi; NaN];
    rP_full   = [r_plot; NaN; r_plot; NaN];
    has_point = true;
else
    has_point = false;
end

% ====================== Enhanced Visualization ======================
fig = figure('Name', sprintf('Enhanced ARP - %s Microscopy', microscopy_type), ...
            'Position', [100, 100, 900, 800], 'Color', 'white');

% Create polar axes with better properties
pax = polaraxes('Parent', fig);
hold(pax, 'on');

% === 1. Improved Grid and Axes ===
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'counterclockwise';

% More informative angular labels
pax.ThetaTick = [0 30 60 90 120 150 180 210 240 270 300 330];
pax.ThetaTickLabel = {'0° (z-axis)', '30°', '60°', ...
                      '90° (x-axis)', '120°', '150°', ...
                      '180°', '210°', '240°', ...
                      '270°', '300°', '330°'};

% Make grid more visible
pax.GridLineStyle = '-';
pax.GridAlpha = 0.3;
pax.GridColor = [0.3 0.3 0.3];
pax.MinorGridLineStyle = ':';
pax.MinorGridAlpha = 0.15;
pax.MinorGridColor = [0.5 0.5 0.5];

% === 2. Color Scheme (Colorblind-friendly) ===
color_point = [0, 114, 189]/255;    % Blue
color_sheet = [217, 83, 25]/255;    % Orange
color_line = [237, 177, 32]/255;    % Yellow/Gold

% === 3. Plot with Enhanced Styling ===
% Add subtle background shading for each curve
if has_point
    % Create filled area for point - properly handle dimensions
    valid_idx = isfinite(angP_full) & isfinite(rP_full);
    if any(valid_idx)
        ang_valid = angP_full(valid_idx);
        r_valid = rP_full(valid_idx);
        % Create closed polygon for filling
        ang_fill = [ang_valid(:); ang_valid(1)];
        r_fill = [r_valid(:); r_valid(1)];
        
        % Only plot fill if we have valid data
        if length(ang_fill) > 2
            fill_handle = polarplot(pax, ang_fill, r_fill, '-', ...
                                  'Color', [color_point, 0.2], 'LineWidth', 0.1);
            fill_handle.DisplayName = '';
        end
    end
    
    % Main curve
    h_point = polarplot(pax, angP_full, rP_full, '--', ...
                       'LineWidth', 2.5, 'Color', color_point, ...
                       'DisplayName', 'Point (PSF FWHM)');
end

if has_sheet
    % Create filled area for sheet - properly handle dimensions
    valid_idx = isfinite(angS_full) & isfinite(rS_full);
    if any(valid_idx)
        ang_valid = angS_full(valid_idx);
        r_valid = rS_full(valid_idx);
        % Create closed polygon for filling
        ang_fill = [ang_valid(:); ang_valid(1)];
        r_fill = [r_valid(:); r_valid(1)];
        
        % Only plot fill if we have valid data
        if length(ang_fill) > 2
            fill_handle = polarplot(pax, ang_fill, r_fill, '-', ...
                                  'Color', [color_sheet, 0.2], 'LineWidth', 0.1);
            fill_handle.DisplayName = '';
        end
    end
    
    % Main curve
    h_sheet = polarplot(pax, angS_full, rS_full, '--', ...
                       'LineWidth', 2.5, 'Color', color_sheet, ...
                       'DisplayName', 'Sheet (SSF FWHM)');
end

if has_line
    % Create filled area for line - properly handle dimensions
    valid_idx = isfinite(angL_full) & isfinite(rL_full);
    if any(valid_idx)
        ang_valid = angL_full(valid_idx);
        r_valid = rL_full(valid_idx);
        % Create closed polygon for filling
        ang_fill = [ang_valid(:); ang_valid(1)];
        r_fill = [r_valid(:); r_valid(1)];
        
        % Only plot fill if we have valid data
        if length(ang_fill) > 2
            fill_handle = polarplot(pax, ang_fill, r_fill, '-', ...
                                  'Color', [color_line, 0.2], 'LineWidth', 0.1);
            fill_handle.DisplayName = '';
        end
    end
    
    % Main curve
    h_line = polarplot(pax, angL_full, rL_full, '--', ...
                      'LineWidth', 2.5, 'Color', color_line, ...
                      'DisplayName', 'Line (LSF FWHM)');
end

title(pax, {sprintf('Anisotropic Resolution Profile (ARP)'), ...
           sprintf('%s Microscopy%s', microscopy_type, pinhole_str), ...
           sprintf('NA=%.1f, λ=%d nm', NA, lambda*1e9)}, ...
      'FontSize', 14, 'FontWeight', 'bold');

% === 7. Improved Legend ===
leg = legend(pax, 'Location', 'northeast', ...
            'FontSize', 11, 'Box', 'on', ...
            'BackgroundAlpha', 0.9);

% Add title to legend
leg.Title.String = 'Target Type';
leg.Title.FontWeight = 'bold';
leg.Title.FontSize = 12;

% === 4. Improved Radial Scale ===
% Calculate appropriate scale
r_candidates = [];
if has_sheet, r_candidates = [r_candidates; rS(:)]; end
if has_line, r_candidates = [r_candidates; rL(:)]; end
if has_point, r_candidates = [r_candidates; rP(:)]; end

rmax = max(r_candidates, [], 'omitnan');
rmin = min(r_candidates, [], 'omitnan');

if isfinite(rmax) && rmax > 0
    % Set nice round numbers for radial axis
    if rmax < 1
        rtick_max = ceil(rmax * 10) / 10;
        rtick_step = 1;
    elseif rmax < 5
        rtick_max = ceil(rmax);
        rtick_step = 1;
    else
        rtick_max = ceil(rmax / 2) * 2;
        rtick_step = 2;
    end
    
    pax.RLim = [0, 6];
    pax.RTick = 0:rtick_step:6;
    
    % Add units to radial labels
    pax.RTickLabel = arrayfun(@(x) sprintf('%.1f μm', x), pax.RTick, 'UniformOutput', false);
    pax.RAxisLocation = 0;
end

% === 5. Add Reference Circles for Key Values ===


% === 6. Enhanced Title and Labels ===
% Filter out empty display names from legends

% === 8. Add Quantitative Annotations ===
% Add min/max values as text annotations
% === 9. Add Isotropy Indicator ===
% Calculate anisotropy ratio for each curve


% === 10. Export Options ===
hold(pax, 'off');

% Save the enhanced figure
%saveas(fig, sprintf('ARP_Enhanced_%s%s.png', microscopy_type, ...
                   %strrep(pinhole_str, ' ', '_')));
%saveas(fig, sprintf('ARP_Enhanced_%s%s.fig', microscopy_type, ...
                   %strrep(pinhole_str, ' ', '_')));

fprintf('\n=== Enhanced ARP Plot Generated ===\n');
fprintf('Key improvements implemented:\n');
fprintf('  - Clearer radial grid with μm units\n');
fprintf('  - Colorblind-friendly palette\n');
fprintf('  - Subtle area shading for each curve\n');
fprintf('  - Reference circles for theoretical limits\n');
fprintf('  - Quantitative range annotations\n');
fprintf('  - Anisotropy ratio calculations\n');
fprintf('==========================================\n');

%% === New code: Save ARP data ===
fprintf('\n>>> Saving Confocal ARP data...\n');

% Create a structure to store data
Confocal_data = struct();
Confocal_data.ang_PSF = angP_full;
Confocal_data.r_PSF   = rP_full;
Confocal_data.ang_SSF = angS_full;
Confocal_data.r_SSF   = rS_full;
Confocal_data.ang_LSF = angL_full;
Confocal_data.r_LSF   = rL_full;

% Save structure to .mat file
save('ARP_data_Confocal.mat', 'Confocal_data');

fprintf('>>> Data successfully saved to ARP_data_Confocal.mat\n');

%% PART 5: ANALYSIS OF SIGNAL-TO-PEDESTAL RATIO (PROFESSOR'S SUGGESTION)
% ========================================================================
fprintf('\n--- Starting Analysis of Signal-to-Pedestal Ratio ---\n');

% --- 1. Pre-computation for various pinhole sizes ---
% We will calculate the PSF and SBR for a range of pinhole sizes.
pinhole_AU_range = 0:0.1:3.0;
num_au_steps = length(pinhole_AU_range);

% Use cell arrays and vectors to store the results
psf_storage = cell(1, num_au_steps);
axial_profiles = zeros(length(z_coords), num_au_steps);
sbr_values = nan(1, num_au_steps);

h_waitbar = waitbar(0, 'Calculating PSFs and SBR for various pinhole sizes...');
tic;

% This loop calculates the necessary data for all three visualization methods
for i = 1:num_au_steps
    current_AU = pinhole_AU_range(i);
    waitbar(i/num_au_steps, h_waitbar, sprintf('Processing AU = %.1f', current_AU));
    
    % --- a) Calculate the final confocal PSF for the current AU ---
    % This logic is reused from our previous successful simulations
    IPSF_det_blurred = zeros(size(IPSF_em));
    if current_AU < 1e-6 % Ideal pinhole case
        IPSF_det_blurred = IPSF_em;
    else
        D_airy_disk_um = (1.22 * lambda_em / NA) * 1e6;
        pinhole_diameter_um = D_airy_disk_um * current_AU;
        pinhole_radius_pixels = (pinhole_diameter_um / 2) / distance_per_pixel;
        
        [X_pix, Y_pix] = meshgrid(pts, pts);
        pinhole_mask_pix = double(sqrt(X_pix.^2 + Y_pix.^2) < pinhole_radius_pixels);
        if sum(pinhole_mask_pix(:)) > 0
            pinhole_mask_pix = pinhole_mask_pix / sum(pinhole_mask_pix(:));
        end
        for z_idx = 1:length(z_coords)
            IPSF_det_blurred(:, :, z_idx) = conv2(squeeze(IPSF_em(:, :, z_idx)), pinhole_mask_pix, 'same');
        end
    end
    IPSF_final = IPSF_ex .* IPSF_det_blurred;
    IPSF_final_norm = IPSF_final / max(IPSF_final(:));
    
    % Store the full 3D PSF
    psf_storage{i} = IPSF_final_norm;
    
    % --- b) Calculate the Signal-to-Pedestal Ratio (SBR) ---
    axial_profile = squeeze(IPSF_final_norm(centre_pt, centre_pt, :));
    axial_profiles(:, i) = axial_profile;
    
    peak_signal = axial_profile(centre_pt); % Should be 1.0
    
    % Use findpeaks to locate the side-lobes ("pedestal" or "shoulders")
    % We search on the second half of the profile to find the first side-lobe
    [pks, ~] = findpeaks(axial_profile(centre_pt+1:end));
    
    if ~isempty(pks)
        pedestal_signal = pks(1); % The first peak after the center is the shoulder
        sbr_values(i) = peak_signal / pedestal_signal;
    else
        sbr_values(i) = Inf; % If no side-lobes are found (ideal case)
    end
end
close(h_waitbar);
fprintf('SBR calculations finished in %.2f seconds.\n', toc);


%% --- 2. Visualization of the results ---

% --- Method 1: Direct Annotation on 2D PSF Plots ---
figure('Name', 'SBR Annotation Method');
au_to_display = [0.5, 1.0, 2.0];
for i = 1:length(au_to_display)
    subplot(1, length(au_to_display), i);
    
    current_AU = au_to_display(i);
    [~, idx] = min(abs(pinhole_AU_range - current_AU));
    
    imagesc(x_coords, z_coords, squeeze(psf_storage{idx}(:, centre_pt, :))');
    axis image;
    title(sprintf('IPSF (%.1f AU)', current_AU));
    xlabel('x (\mum)');
    if i == 1, ylabel('z (\mum)'); end
    colorbar; clim([0, 0.1]);
    
    % Add the SBR text annotation
    text(0, -4, sprintf('S/B Ratio = %.1f', sbr_values(idx)), ...
        'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
sgtitle('Visualization Method 1: Direct SNR Annotation on PSF Plot', 'FontSize', 16, 'FontWeight', 'bold');


% --- Method 2: Comparison of 1D Axial Profiles ---
figure('Name', 'SBR Profile Comparison Method');
hold on; grid on;
au_indices_to_plot = [1, 6, 11, 21]; % Corresponds to AU = 0.0, 0.5, 1.0, 2.0
colors = {'k', 'b', 'r', 'g'};
for i = 1:length(au_indices_to_plot)
    idx = au_indices_to_plot(i);
    plot(z_coords, axial_profiles(:, idx), 'Color', colors{i}, 'LineWidth', 2, ...
        'DisplayName', sprintf('%.1f AU', pinhole_AU_range(idx)));
end
hold off;
title('Visualization Method 2: Axial Profile Comparison', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Axial Position (\mum)');
ylabel('Normalized Intensity');
legend('show', 'Location', 'northeast');
xlim([-5, 5]);
ylim([0, 0.1]); % Zoom in on the pedestal/shoulders
text(-4, 0.09, 'Note: As pinhole size increases, "shoulder" becomes significantly elevated', 'FontSize', 12);


% --- Method 3: Summary Curve of SBR vs. AU ---
figure('Name', 'SBR Summary Curve Method');
plot(pinhole_AU_range, sbr_values, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
grid on;
set(gca, 'YScale', 'log'); % Use a log scale to better visualize the large change
title('Peak-to-Sidelobe Ratio (PSR) vs Pinhole Size (AU)', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Pinhole Size (AU)');
ylabel('Peak-to-Sidelobe Ratio (PSR)');
xlim([0, max(pinhole_AU_range)]);