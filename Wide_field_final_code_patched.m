clear; close all; clc;
%% Basic settings
NA = 0.5; % NA of lens
n = 1.33; % refractive index
lambda = 500e-9; % wavelength in metres%
wavenumber = n /lambda; % wavenumber in metres^-1
%wavenumber = 2*pi/lambda;

%--------------
grid_size = 128; % number of points in grid in 1D
pts = [-grid_size:1:grid_size-1]; % 1D array of grid points
centre_pt = grid_size + 1; % number of the element corresponding to the centre of the array
[X, Y, Z] = meshgrid(pts, pts, pts);
R = sqrt(X.^2+Y.^2+Z.^2); % Circle Ball

%---------------
shell_radius = grid_size/2;
shell_thickness = 1; % thickness of shell in voxels
shell = shell_radius-shell_thickness/2 < R & R <= shell_radius+shell_thickness/2;

spatial_frequency_per_pixel_um = (wavenumber/shell_radius) * 1e-6;
spatial_frequency_range = -grid_size*spatial_frequency_per_pixel_um:spatial_frequency_per_pixel_um:(grid_size-1)*spatial_frequency_per_pixel_um;

X_scaled = X*spatial_frequency_per_pixel_um;
Y_scaled = Y*spatial_frequency_per_pixel_um;
Z_scaled = Z*spatial_frequency_per_pixel_um;

fc = (2*NA/lambda) * 1e-6;
fprintf('fc = %f um\n', fc);

phi = atan2(sqrt(X.^2+Y.^2), Z);
cone = (0 <= phi & phi <= asin(NA/n));
CTF = double(cone & shell); % generate the cap by AND of shell and cone
APSF = ifftshift(ifftn(fftshift(CTF))); % calculate the amplitude PSF
IPSF = abs(APSF).^2; % calculate the intensity PSF
OTF = ifftshift(fftn(fftshift(IPSF))); % calculate the OTF
MTF = abs(OTF); % calculate the MTF
MTF = MTF/max(MTF(:)); % normalise the MTF
IPSF = IPSF./max(IPSF(:));% normalise the IPSF

% Unit conversion between frequency and spatial domains (generate real-space scale in meters and micrometers)
distance_per_pixel = 1/(spatial_frequency_per_pixel_um*2*grid_size);
distance_range = -grid_size*distance_per_pixel:distance_per_pixel:(grid_size-1)*distance_per_pixel;
dist_um = distance_range ;
[X_real, Y_real, Z_real] = meshgrid(dist_um, dist_um, dist_um);
x_coords = dist_um; y_coords = dist_um; z_coords = dist_um;
xx = 0.51*lambda/NA * 1e6;
zz = 0.89 * lambda / (n - sqrt(n^2 - NA^2))* 1e6;
spatial_freq_coords = spatial_frequency_range;



%% Figures: 3D、2D visualization of the CTF (pupil-cone)

% Figure 1: 3D visualization of the CTF (pupil-cone)
figure
CTF_iso_surface = isosurface(X_scaled, Y_scaled, Z_scaled, CTF, 0.5);
p1 = patch(CTF_iso_surface);
isonormals(X_scaled, Y_scaled, Z_scaled, CTF, p1)

p1.FaceColor = 'blue';
p1.EdgeColor = 'none';
daspect([1,1,1])
view(3), axis equal
camlight; lighting gouraud
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('3D Plot of Pupil Collection Cone (CTF)', 'FontSize', 24, 'FontName','Times New Roman', 'FontWeight', 'bold');
xlabel('u_x (cycles/um)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
ylabel('u_y (cycles/um)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
zlabel('u_z (cycles/um)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');


%zlim([-fc, fc])

%legend('CTF (arbitrary scaling on x,y,z axes)')

% Figure 2: 2D view of the CTF (u_y=0, x-z plane only)
figure
imagesc(spatial_frequency_range,spatial_frequency_range, squeeze(CTF(centre_pt, :, :))')
axis xy
axis equal
axis tight
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('CTF, u_y = 0', 'FontSize', 24, 'FontName', 'Times New Roman','FontWeight', 'bold');
ylabel('u_z (cycles/um)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlabel('u_x (cycles/um)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
%ylim([-1.1*fc, 1.1*fc])
colorbar

figure
imagesc(spatial_frequency_range,spatial_frequency_range, squeeze(CTF(:, :, centre_pt+60)))
axis equal
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('CTF, u_z = 0', 'FontSize', 24, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlabel('u_x (cycles/um)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
ylabel('u_y (cycles/um)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
clim([0, 1])
xlim([-5, 5])
colorbar

%% Figures: 3D、2D visualization of the MTF

% Figure 3: 3D isosurface visualization of the MTF
figure
MTF_iso_surface = isosurface(X_scaled, Y_scaled, Z_scaled, MTF, 0.001);
p = patch(MTF_iso_surface);
isonormals(X_scaled, Y_scaled, Z_scaled, MTF, p);

p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3), axis equal
%view(52, 28), axis equal
camlight; lighting gouraud;
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('MTF Boundary Isosurface', 'FontSize', 24, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlabel('u_x (cycles/µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
ylabel('u_y (cycles/µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
zlabel('u_z (cycles/µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
%zlim([-2*fc, 2*fc])

title('MTF boundary isosurface')

% Figure 4,5: 2D MTF slice (u_y=0)  "x-z plane"
figure
imagesc(spatial_frequency_range, spatial_frequency_range, squeeze(MTF(centre_pt, :, :))')
axis equal
axis tight
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('MTF, u_y = 0, range =[0,1]', 'FontSize', 24, 'FontName', 'Times New Roman');
ylabel('u_z (cycles/µm)', 'FontSize', 20, 'FontName', 'Times New Roman');
xlabel('u_x (cycles/µm)', 'FontSize', 20, 'FontName', 'Times New Roman');
xlim([-2*fc, 2*fc]); ylim([-2*fc, 2*fc]);
set(gca, 'FontName', 'Times New Roman');
%xlim([-2*fc, 2*fc])
%ylim([-2*fc, 2*fc])
colorbar
clim([0, 1])

% MTF main analysis figure
figure
imagesc(spatial_frequency_range, spatial_frequency_range, squeeze(MTF(centre_pt, :, :))')
axis equal
axis tight
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('MTF, u_y = 0, range = [0, 0.1]', 'FontSize', 24, 'FontName', 'Times New Roman','FontWeight', 'bold');
ylabel('u_z (cycles/µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlabel('u_x (cycles/µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlim([-2*fc, 2*fc]); ylim([-1.5*fc, 1.5*fc]);
colorbar
clim([0, 0.1])

hold on; % Allow drawing on top of the existing image

% Get the current y-axis limits to draw the lines across the full height
yLims = ylim;

% Draw vertical dashed lines at the positive and negative cutoff frequencies
plot([fc, fc], yLims, 'r--', 'LineWidth', 1, 'DisplayName', 'Cutoff Frequency');
plot([-fc, -fc], yLims, 'r--', 'LineWidth', 1, 'HandleVisibility', 'off'); % Hide from legend

% Add text labels next to the lines for clarity
text(fc + 0.1, 0.5, sprintf('f_c = %.2f um', fc), ...
    'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');
text(-fc - 1.5, 0.5, sprintf('-f_c = -%.2f um', fc), ...
    'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');

hold off; % Finish drawing on the plot

% Other MTF figure

% Figure 6,7: 2D MTF slice (u_z=0)  "x-y plane"
figure
imagesc(spatial_frequency_range, spatial_frequency_range, squeeze(MTF(:, :, centre_pt)))
axis equal
axis tight
title('MTF, u_z = 0, range = [0, 1]', 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('u_y (cycles/µm)', 'FontSize', 16, 'FontName', 'Times New Roman');
ylabel('u_x (cycles/µm)', 'FontSize', 16, 'FontName', 'Times New Roman');
xlim([-2*fc, 2*fc]); ylim([-2*fc, 2*fc]);
set(gca, 'FontName', 'Times New Roman');
%xlim([-2*fc, 2*fc])
%ylim([-2*fc, 2*fc])
clim([0, 1])
colorbar

figure
imagesc(spatial_frequency_range, spatial_frequency_range, squeeze(MTF(:, :, centre_pt)))
axis equal
axis tight
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('MTF, u_z = 0, range = [0, 0.1]', 'FontSize', 24, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlabel('u_y (cycles/µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
ylabel('u_x (cycles/µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlim([-2*fc, 2*fc]); ylim([-2*fc, 2*fc]);
%xlim([-2*fc, 2*fc])
%ylim([-2*fc, 2*fc])
clim([0, 0.1])
colorbar


%% Figures: 3D、2D visualization of the IPSF, Lateral and Axial resolution calculation

% Figure 8: 3D_IPSF

figure
IPSF_iso_surface = isosurface(X_real, Y_real, Z_real, IPSF, 0.001);
p2 = patch(IPSF_iso_surface);
isonormals(X_real, Y_real, Z_real, IPSF, p2)

p2.FaceColor = 'blue';
p2.EdgeColor = 'none';
daspect([1,1,1])
view(3), axis equal
camlight; lighting gouraud
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('3D Plot of IPSF', 'FontSize', 24, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlabel('x (µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
ylabel('y (µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
zlabel('z (µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');



% FIgure 9,10,11: 2D IPSF slice

xx = 0.51 * lambda / NA * 1e6; 
zz = 0.89 * lambda /(n - sqrt(n^2 - NA^2))* 1e6; 

%── (1) (u_z=0)  "x-y plane" ──
figure
IPSF_xy = squeeze(IPSF(:, :, centre_pt));  % 大小 256×256
imagesc(dist_um, dist_um, IPSF_xy');
axis image
xlim([-1.5*xx, 1.5*xx])
ylim([-1.5*xx, 1.5*xx])
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('IPSF Slice at z = 0 (x-y plane)', 'FontSize', 24, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlabel('y (µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
ylabel('x (µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
colorbar

% ── (2) (u_y=0)  "x-z plane" ──
figure
IPSF_xz = squeeze(IPSF(centre_pt, :, :));  
imagesc(dist_um, dist_um, IPSF_xz');
axis image
xlim([-1.5*zz, 1.5*zz])
ylim([-1.5*zz, 1.5*zz])
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('IPSF Slice at y = 0 (x-z plane)', 'FontSize', 24, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlabel('x (µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
ylabel('z (µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
colorbar


%── (3) (u_x=0)  "y-z plane" ──
figure
IPSF_yz = squeeze(IPSF(:, centre_pt, :));  
imagesc(dist_um, dist_um, IPSF_yz');
axis image
xlim([-1.5*zz, 1.5*zz])
ylim([-1.5*zz, 1.5*zz])
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



% Figure 12,13: Lateral profile and FWHM (version with FWHM value annotation)
xx = 0.51*lambda/NA * 1e6; 
figure;
lateral_line = squeeze( IPSF(centre_pt,:,centre_pt) );
x_coords      = dist_um; 
factor = 20000;
xp = linspace( x_coords(1), x_coords(end), numel(x_coords)*factor );
yp = interp1( x_coords, lateral_line./max(lateral_line), xp, 'spline' );
plot(xp, yp, 'b-', 'LineWidth',2, 'DisplayName', 'PSF Profile');
hold on;
grid on;
xlim([-2*zz, 2*zz]);
ylim([0, 1.05]);

% --- Calculate FWHM ---
halfVal = 0.5; 
[~, peakIdx1] = max( yp );
leftIdx1 = find( yp(1:peakIdx1) < halfVal, 1, 'last' );
rightIdx1 = find( yp(peakIdx1:end) < halfVal, 1, 'first' ) + peakIdx1 - 1;
xL1 = interp1( yp([leftIdx1,   leftIdx1+1]), xp([leftIdx1,   leftIdx1+1]), halfVal );
xR1 = interp1( yp([rightIdx1-1, rightIdx1]), xp([rightIdx1-1, rightIdx1]), halfVal );
FWHM_lateral_um = abs(xR1 - xL1);

% --- Draw FWHM indicator lines ---
plot([xL1, xR1], [0.5, 0.5], 'r--', 'LineWidth', 2, 'DisplayName', 'FWHM');

% [NEW] Annotate FWHM value on the plot
text_x_pos = (xL1 + xR1) / 2; % x coordinate at center of dashed line
text_y_pos = 0.4;             % y coordinate below the dashed line
fwhm_string = sprintf('FWHM = %.3f µm', FWHM_lateral_um);
text(text_x_pos, text_y_pos, fwhm_string, ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 18, ...
    'FontName', 'Times New Roman',...
    'FontWeight', 'bold',...
    'Color', 'r');

% --- Format title and legend ---
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('Lateral PSF Profile', 'FontSize', 24, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlabel('Lateral Position (µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold'); 
ylabel('Normalized Intensity', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
legend('show', 'Location', 'Best', 'FontName', 'Times New Roman', 'FontSize', 16);
hold off;

fprintf('Theory Lateral FWHM after spline interpolation = %.4f µm\n', xx);
fprintf('Measure Lateral FWHM after spline interpolation = %.4f µm\n', FWHM_lateral_um);

% Figure 14,15: Axial profile and FWHM

zz = 0.89 * lambda / (n - sqrt(n^2 - NA^2))* 1e6; 
figure;
axial_line = squeeze( IPSF( centre_pt, centre_pt, : ) );  
z_coords   = dist_um;
factor = 20;
xp2 = linspace( z_coords(1), z_coords(end), numel(z_coords)*factor );
yp2 = interp1( z_coords, axial_line./max(axial_line), xp2, 'spline' );
plot(xp2, yp2, 'b-', 'LineWidth',2, 'DisplayName', 'PSF Profile');
hold on;
grid on;
xlim([-2*zz, 2*zz]);
ylim([0, 1.05]);

% --- Calculate FWHM ---
halfVal = 0.5;
[~, peakIdx2] = max( yp2 );
leftIdx2 = find( yp2(1:peakIdx2) < halfVal, 1, 'last' );
rightIdx2 = find( yp2(peakIdx2:end) < halfVal, 1, 'first' ) + peakIdx2 - 1;
xL2 = interp1( yp2([leftIdx2,   leftIdx2+1]), xp2([leftIdx2,   leftIdx2+1]), halfVal );
xR2 = interp1( yp2([rightIdx2-1, rightIdx2]), xp2([rightIdx2-1, rightIdx2]), halfVal );
FWHM_axial_um = abs(xR2 - xL2);

% --- Draw FWHM indicator lines ---
plot([xL2, xR2], [0.5, 0.5], 'r--', 'LineWidth', 2, 'DisplayName', 'FWHM');

% [NEW] Annotate FWHM value on the plot
text_x_pos = (xL2 + xR2) / 2; % x coordinate at center of dashed line
text_y_pos = 0.4;             % y coordinate below the dashed line
fwhm_string = sprintf('FWHM = %.3f µm', FWHM_axial_um);
text(text_x_pos, text_y_pos, fwhm_string, ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 18, ...
    'FontName', 'Times New Roman',...
    'FontWeight', 'bold',...
    'Color', 'r');

% --- Format title and legend ---
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('Axial PSF Profile', 'FontSize', 24, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlabel('Axial Position (µm)', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold'); 
ylabel('Normalized Intensity', 'FontSize', 20, 'FontName', 'Times New Roman','FontWeight', 'bold');
legend('show', 'Location', 'Best', 'FontName', 'Times New Roman', 'FontSize', 16);
hold off;

fprintf('Theory Axial FWHM after spline interpolation = %.4f µm\n', zz);
fprintf('Measure Axial FWHM after spline interpolation = %.4f µm\n', FWHM_axial_um);


%% ==== FWHM fraction (3D) + optional axial FWHM ====
% Assumption: IPSF is already a 3D intensity kernel (can be normalized or unnormalized - this metric is amplitude-scaling insensitive)
%       dist_um / X_real, Y_real, Z_real are already defined (real space coordinates)

% --- Voxel size and volume (µm) ---
dx = dist_um(2) - dist_um(1);
dy = dx; dz = dx;                   % your grid is equally spaced
voxelVol = dx*dy*dz;                % voxel volume (µm^3)

% --- 3D FWHM volume energy fraction ---
V = IPSF;                           % intensity volume
Imax = max(V(:));
mask3D = (V >= 0.5*Imax);           % FWHM mask
E_total = sum(V(:)) * voxelVol;     % total energy
E_fwhm  = sum(V(mask3D)) * voxelVol;% energy within FWHM
frac_fwhm3D = E_fwhm / E_total;     % <-- the metric you want

% --- (Optional) 1D FWHM in x/y/z directions (linear interpolation) ---
% Extract three profiles through the voxel center
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
IPSF_iso_surface = isosurface(X_real, Y_real, Z_real, IPSF, 0.01);
p2 = patch(IPSF_iso_surface);
isonormals(X_real, Y_real, Z_real, IPSF, p2)
p2.FaceColor = 'blue';
p2.EdgeColor = 'none';
daspect([1,1,1])
view(3), axis equal
camlight; lighting gouraud
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('3D Plot of Wide-field Microscope IPSF', 'FontSize', 24, 'FontName', 'Times New Roman','FontWeight', 'bold');
xlabel('x (µm)', 'FontSize', 16, 'FontName', 'Times New Roman','FontWeight', 'bold');
ylabel('y (µm)', 'FontSize', 16, 'FontName', 'Times New Roman','FontWeight', 'bold');
zlabel('z (µm)', 'FontSize', 16, 'FontName', 'Times New Roman','FontWeight', 'bold');


% ===== Add 3D FWHM fraction information in the upper right corner =====
txt_str = sprintf('3D FWHM fraction = %.3f\nFWHM x/y/z = %.2f/%.2f/%.2f µm', ...
                  frac_fwhm3D, FWHM_x, FWHM_y, FWHM_z);

% Calculate upper right corner position (slightly outside data range)
x_pos = max(x_coords)*0.9;
y_pos = max(y_coords)*0.9;
z_pos = max(z_coords)*(-0.25);

t = text(x_pos, y_pos, z_pos, txt_str, ...
         'FontSize', 16, 'FontWeight', 'bold', ...
         'FontName', 'Times New Roman', ...
         'BackgroundColor', [1 1 1 0.75], ... % 白色+半透明
         'EdgeColor', [0 0 0], ...
         'Margin', 5, ...
         'VerticalAlignment', 'top', ...
         'HorizontalAlignment', 'left');


function w = local_fwhm_1d(x, y)
    % 线性插值的一维 FWHM（y 可未归一化）
    x = x(:); y = y(:);
    [ymax, imax] = max(y);
    half = ymax/2;

    % 左半最大到半高的交点
    iL = find(y(1:imax) <= half, 1, 'last');
    if isempty(iL) || iL==imax
        xL = x(1);
    else
        xL = interp1(y([iL, iL+1]), x([iL, iL+1]), half, 'linear');
    end

    % 右半从峰到半高的交点
    iRrel = find(y(imax:end) <= half, 1, 'first');
    if isempty(iRrel) || imax-1+iRrel==imax
        xR = x(end);
    else
        iR = imax-1+iRrel;
        xR = interp1(y([iR-1, iR]), x([iR-1, iR]), half, 'linear');
    end

    w = abs(xR - xL);
end

%% ===== Component 1: Quantitative metrics (FWHM_x/y/z, HMEF, Tail energy, radial MTF) =====
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
    [fc_thresh, fc_mtf10, fc_mtf50] = mtf_cutoff_linear(fr_centers(:), mtf_rad(:), 1e-3);
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


%% Theory Research 1: NA vs. theory lateral and axial resolution (T.Wilson's Paper)

% Batch scan NA. From 0.5 to 1.3, with a step size set to 0.1.
NA_list = 0.2 : 0.1 : 1;
nNA     = numel(NA_list);

% Pre allocation
meas_lat       = nan(nNA,1);
meas_ax        = nan(nNA,1);
theo_lat       = nan(nNA,1);
theo_ax_ray    = nan(nNA,1);
theo_ax_non    = nan(nNA,1);

n = 1; % refractive index
lambda = 500e-9; % wavelength in metres%
wavenumber = n /lambda; % wavenumber in metres^-1
%wavenumber = 2*pi/lambda;

%--------------
grid_size = 256; % number of points in grid in 1D
pts = [-grid_size:1:grid_size-1]; % 1D array of grid points
centre_pt = grid_size + 1; % number of the element corresponding to the centre of the array
[X, Y, Z] = meshgrid(pts, pts, pts);
R = sqrt(X.^2+Y.^2+Z.^2); % Circle Ball

%---------------
shell_radius = grid_size/2;
shell_thickness = 1; % thickness of shell in voxels
shell = shell_radius-shell_thickness/2 < R & R <= shell_radius+shell_thickness/2;

spatial_frequency_per_pixel_um = (wavenumber/shell_radius) * 1e-6;
spatial_frequency_range = -grid_size*spatial_frequency_per_pixel_um:spatial_frequency_per_pixel_um:(grid_size-1)*spatial_frequency_per_pixel_um;

X_scaled = X*spatial_frequency_per_pixel_um;
Y_scaled = Y*spatial_frequency_per_pixel_um;
Z_scaled = Z*spatial_frequency_per_pixel_um;

for i = 1 : nNA
    NA = NA_list(i);
    phi = atan2(sqrt(X.^2+Y.^2), Z);
    cone = (0 <= phi & phi <= asin(NA/n));
    CTF = double(cone & shell); 

    APSF = ifftshift(ifftn(fftshift(CTF))); % calculate the amplitude PSF
    IPSF = abs(APSF).^2; % calculate the intensity PSF
    OTF = ifftshift(fftn(fftshift(IPSF))); % calculate the OTF
    MTF = abs(OTF); % calculate the MTF
    MTF = MTF/max(MTF(:)); % normalise the MTF
    IPSF = IPSF./max(IPSF(:));% normalise the IPSF

    distance_per_pixel = 1/(spatial_frequency_per_pixel_um*2*grid_size);
    distance_range = -grid_size*distance_per_pixel:distance_per_pixel:(grid_size-1)*distance_per_pixel;
    dist_um = distance_range ;
    [X_real, Y_real, Z_real] = meshgrid(dist_um, dist_um, dist_um);

    lateral_line = squeeze( IPSF(centre_pt,:,centre_pt) );
    x_coords      = dist_um; 
    factor = 20000;  % 
    xp = linspace( x_coords(1), x_coords(end), 256*factor );
    yp = interp1( x_coords, lateral_line./max(lateral_line), xp, 'spline' );

    xp_lat{i}=xp; yp_lat{i}=yp;

    halfVal = 0.5; 
    [~, peakIdx1] = max( yp );
    leftIdx1 = find( yp(1:peakIdx1) < halfVal, 1, 'last' );
    rightIdx1 = find( yp(peakIdx1:end) < halfVal, 1, 'first' ) + peakIdx1 - 1;

    xL1 = interp1( yp([leftIdx1,   leftIdx1+1]), ...
              xp([leftIdx1,   leftIdx1+1]), ...
              halfVal );

    xR1 = interp1( yp([rightIdx1-1, rightIdx1]), ...
              xp([rightIdx1-1, rightIdx1]), ...
              halfVal );
    meas_lat(i) = abs(xR1 - xL1);

%%-----%%
    axial_line = squeeze( IPSF( centre_pt, centre_pt, : ) );  
    z_coords   = dist_um; 
    factor = 20;  % 
    xp2 = linspace( z_coords(1), z_coords(end), 256*factor );
    yp2 = interp1( z_coords, axial_line./max(axial_line), xp2, 'spline' );

    halfVal = 0.5; 
    [~, peakIdx2] = max( yp2 );
    leftIdx2 = find( yp2(1:peakIdx2) < halfVal, 1, 'last' );
    rightIdx2 = find( yp2(peakIdx2:end) < halfVal, 1, 'first' ) + peakIdx2 - 1;

    xL2 = interp1( yp2([leftIdx2,   leftIdx2+1]), ...
              xp2([leftIdx2,   leftIdx2+1]), ...
              halfVal );

    xR2 = interp1( yp2([rightIdx2-1, rightIdx2]), ...
              xp2([rightIdx2-1, rightIdx2]), ...
              halfVal );
    meas_ax(i) = abs(xR2 - xL2);

    %%--- %%
    theo_lat(i)    = 0.51 * lambda/NA * 1e6;                      % Airy Lateral Resolution
    theo_ax_ray(i) = 1.77 * n*lambda/(NA^2) * 1e6;                % Rayleigh Axial Resolution
    theo_ax_non(i) = 0.89 * lambda/(n - sqrt(n^2 - NA^2)) * 1e6;  % Axial Resolution from T.Wilson's paper
end

T = table( NA_list(:), meas_lat, theo_lat, meas_ax, theo_ax_ray, theo_ax_non,...
    'VariableNames', { ...
      'NA', ...
      'Measured_Lateral_um', 'Theo_Lateral_um',  'Measured_Axial_um', ...
      'Theo_Axial_Ray_um', 'Theo_Axial_Nonparaxial_um' } );
disp(T);
writetable(T, 'NA_scan_results.csv');


%% Plot Part


%--- Plot 1: Lateral FWHM vs NA ---
figure;
hold on;
plot(NA_list, meas_lat,'b-o','LineWidth',2, 'DisplayName', 'Measured');
plot(NA_list, theo_lat,'r--','LineWidth',2, 'DisplayName', 'Theory');
hold off;
grid on;
% --- Formatting ---
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold','FontSize', 14);
title('Lateral FWHM vs. NA', 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Numerical Aperture (NA)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Lateral FWHM (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
lgd = legend('show', 'Location', 'Best');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');



%--- Plot 2: Axial FWHM vs NA ---
figure;
hold on;
plot(NA_list, meas_ax,'b-o','LineWidth',2, 'DisplayName', 'Measured');
plot(NA_list, theo_ax_ray,'r--','LineWidth',2, 'DisplayName', 'Common Theory');
plot(NA_list, theo_ax_non,'g--','LineWidth',2, 'DisplayName', 'T.Wilson Theory');
hold off;
grid on;
% --- Formatting ---
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold','FontSize', 14);
title('Axial FWHM vs. NA', 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Numerical Aperture (NA)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Axial FWHM (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
lgd = legend('show', 'Location', 'Best');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');
xlim([0.4, 1]);

%--- [NEW] Plot 3: Detailed Lateral Profile Comparison (for last NA) ---
% Use data from the last iteration of the loop for detailed plotting
figure;
hold on;
grid on;
% 1. Draw PSF profile curve
plot(xp, yp, 'b-', 'LineWidth', 2, 'DisplayName', 'PSF lateral profile (interped)');
% 2. Draw y=0.5 half-maximum line
plot(xlim, [0.5 0.5], 'k--', 'LineWidth', 2, 'DisplayName', 'Half peak line y = 0.5');
% 3. Draw theoretical FWHM (green dashed line)
th_lat_half = theo_lat(end)/2;
plot([-th_lat_half, -th_lat_half], [0, 0.5], 'g--', 'LineWidth', 2, 'DisplayName', 'Theoretical FWHM (green)');
plot([th_lat_half, th_lat_half], [0, 0.5], 'g--', 'LineWidth', 2, 'HandleVisibility', 'off');
% 4. Draw measured FWHM (red solid line)
plot([xL1, xL1], [0, 0.5], 'r-', 'LineWidth', 2, 'DisplayName', 'Measured FWHM (red)');
plot([xR1, xR1], [0, 0.5], 'r-', 'LineWidth', 2, 'HandleVisibility', 'off');
% 5. Draw red dots at half-maximum points
plot([xL1, xR1], [0.5, 0.5], 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'HandleVisibility', 'off');
hold off;
% --- Formatting ---
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold','FontSize', 14);
title(sprintf('IPSF Lateral Profile (NA = %.1f)', NA_list(end)), 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('x (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Normalized Intensity', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
lgd = legend('show', 'Location', 'Best');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');
xlim([-1, 1]);
ylim([0, 1.05]);


%--- [NEW] Plot 4: Detailed Axial Profile Comparison (for last NA) ---
% Use data from the last iteration of the loop for detailed plotting
figure;
hold on;
grid on;
% 1. Draw PSF profile curve
plot(xp2, yp2, 'b-', 'LineWidth', 2, 'DisplayName', 'PSF axial profile (interped)');
% 2. Draw y=0.5 half-maximum line
plot(xlim, [0.5 0.5], 'k--', 'LineWidth', 2, 'DisplayName', 'Half peak line y = 0.5');
% 3. Draw theoretical FWHM (green dashed line)
th_ax_half = theo_ax_non(end)/2;
plot([-th_ax_half, -th_ax_half], [0, 0.5], 'g--', 'LineWidth', 2, 'DisplayName', 'Theoretical FWHM (green)');
plot([th_ax_half, th_ax_half], [0, 0.5], 'g--', 'LineWidth', 2, 'HandleVisibility', 'off');
% 4. Draw measured FWHM (red solid line)
plot([xL2, xL2], [0, 0.5], 'r-', 'LineWidth', 2, 'DisplayName', 'Measured FWHM (red)');
plot([xR2, xR2], [0, 0.5], 'r-', 'LineWidth', 2, 'HandleVisibility', 'off');
% 5. Draw red dots at half-maximum points
plot([xL2, xR2], [0.5, 0.5], 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'HandleVisibility', 'off');
hold off;
% --- Formatting ---
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold','FontSize', 14);
title(sprintf('IPSF Axial Profile (NA = %.1f)', NA_list(end)), 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('z (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Normalized Intensity', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
lgd = legend('show', 'Location', 'Best');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlim([-1, 1]);

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
%% Sheet Object(HORIZONTAL) Figures

% Create a sheet object and visualize
% Create a 3D sheet: one-voxel thick at z=0, uniform in x/y
sheet3D_z0 = zeros(size(IPSF));      % size = [Nx,Ny,Nz]
[~, iz0] = min(abs(dist_um));     % find index of z=0
sheet3D_z0(:,:,iz0) = 1;             % set that plane to 1
% Visualize the z=0 sheet
figure;
p = patch( isosurface(X_real, Y_real, Z_real, sheet3D_z0, 0.5) );
isonormals(X_real, Y_real, Z_real, sheet3D_z0, p);
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
I_freq3D = I_freq3D./max(I_freq3D(:));% normalise the IPSF

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
imagesc(x_coords, z_coords, squeeze(real(I_freq3D(centre_pt,:,:)))');
axis image; 
set(gca,'FontName', 'Times New Roman','FontSize', 14);
title('b) Image after PSF convolution');
xlabel('x (\mum)');
ylabel('z (\mum)')
colorbar

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
% --- Draw frequency domain view of horizontal sheet (Red: Sheet Spectrum, Blue: OTF) ---
figure;

FT_sheet_z0_centered = fftshift(FT_sheet);
otf_slice_2d        = squeeze(MTF(centre_pt, :, :));                % [uy,ux,uz] -> fix uy
sheet_spec_2d       = abs( squeeze( FT_sheet_z0_centered(centre_pt, :, :) ) );
log_otf     = log(1 + otf_slice_2d);
log_otf_norm = log_otf / max(log_otf(:));
sheet_norm = sheet_spec_2d / max(sheet_spec_2d(:));  % normalize to [0,1]
mask_line  = sheet_norm > 0.2;                       % threshold can be adjusted as needed

imagesc(spatial_freq_coords, spatial_freq_coords, log_otf_norm');
axis xy image;
title('Frequency View for Horizontal Sheet', 'FontSize', 16, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Axial Spatial Frequency u_z (\mum^{-1})', 'FontSize', 16, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Lateral Spatial Frequency u_x (\mum^{-1})', 'FontSize', 16, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
lgd = legend('show', 'Location', 'Best');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold', 'TextColor', 'white');
colormap(parula);          % deep blue to yellow
clim([0,0.01]);           % MTF maximum value approximately 0.01, fine-tune as needed
hcb = colorbar;
hcb.Label.String = 'MTF value';hold on;
[rows, cols] = find(mask_line);
u_z_line = spatial_freq_coords(rows);  % rows correspond to y axis
u_x_line = spatial_freq_coords(cols);  % columns correspond to x axis
plot(u_z_line, u_x_line, '.r', 'MarkerSize', 8);
xlim([-2*fc, 2*fc]);
ylim([-2*fc, 2*fc]);
set(gca,'Layer','top');
hold off;

%% === Sheet metrics figure: Axial response with FWHM_z & TER ===
% Calculate z-direction energy distribution Ez after volume energy normalization (directly based on I_freq3D)
I_abs = abs(I_freq3D);
E_total = sum(I_abs(:)) + eps;
Ez = squeeze(sum(sum(I_abs,1),2)) / E_total;   % energy per z-plane (sum = 1)

z0_um = 2 * FWHM_f;                             % commonly used threshold in literature: |z|>=2·FWHM_z
tail_mask = (abs(z_coords) >= z0_um);
TER_sheet = sum(Ez(tail_mask));                 % tail energy ratio (smaller is better)

fprintf('Sheet metrics: FWHM_z = %.3f um | TER(|z|>=%.3f um) = %.4f\n', ...
        FWHM_f, z0_um, TER_sheet);

% Create a separate "metrics version" response plot (clearer)
figure('Name','Sheet Axial Response (metrics)','Color','w');
plot(z_coords, Ez/max(Ez), 'b-', 'LineWidth', 1.8); hold on; grid on;
yline(0.5,':','Half-max');
% Mark FWHM endpoints on z axis (convert your zp/ifq FWHM results to nearest grid points)
[~,iL] = min(abs(z_coords - (xL_f))); 
[~,iR] = min(abs(z_coords - (xR_f)));
xline(z_coords(iL),':'); xline(z_coords(iR),':');
title(sprintf('Sheet axial response: FWHM_z=%.3f um, TER_{|z|>=%.3f um}=%.3f', ...
      FWHM_f, z0_um, TER_sheet));
xlabel('z (um)'); ylabel('Normalized energy/profile');
legend('z-profile (energy projection)','Location','best');



%% Sheet Object(ROTATE) Figures

% 1) Construct arbitrary angle thin sheet (rotate around Y axis) and visualize

theta = 0; % define rotation angle (60 degrees)
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
% 1) FWHM_n: directly based on the normal direction profile you already obtained (line_coords_1d, profile_freq_norm)
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

% 2) TER_normal: integrate 3D intensity over planes parallel to the sheet to get normal direction energy distribution E(u)
I_abs_s = abs(I_freq_slanted);                    % 数值稳健
E_tot_s = sum(I_abs_s(:)) + eps;

% Construct a set of basis vectors orthogonal to normal_vec: t1_vec, t2_vec (take any direction not collinear with normal and orthogonalize)
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

% 3) Metrics plot (along normal direction)
figure('Name','Slanted sheet axial response (along normal)','Color','w');
plot(u_coords, E_u / max(E_u), 'b-', 'LineWidth', 1.8); hold on; grid on;
yline(0.5,':','Half-max');
% Mark FWHM endpoints on normal coordinates
if isfinite(FWHM_n)
    xline(-FWHM_n/2,':'); xline(FWHM_n/2,':');
end
title(sprintf('Along-normal response: FWHM_n=%.3f um, TER_{|u|>=%.3f}=%.3f (\\theta=%.1f^\\circ)', ...
      FWHM_n, u0, TER_normal, theta*180/pi));
xlabel('Distance along normal u (um)'); ylabel('Normalized energy/profile');
xlim([-10,10]);
ylim([0,1]);

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
clim([0.9 1]);

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

% Frequency spectrum visualization (MODIFIED TO DRAW A SMOOTH ANALYTICAL LINE)

figure;
% 1) Draw OTF background (exactly the same as before)
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
u_z_line_smooth = t * line_direction(2); % z 坐标
u_x_line_smooth = t * line_direction(1); % x 坐标
plot(u_x_line_smooth, u_z_line_smooth, 'r--', 'LineWidth', 2); % Use '--' to draw dashed line
xlim([-2*fc, 2*fc]);
ylim([-2*fc, 2*fc]);
hold off;







%% Sheet Object Resolution Figure (Any Angles -> FWHM)

fprintf('\nStarting analysis of FWHM vs. Sheet Angle...\n');

% --- 1. Initialize parameters ---
N_angles = 91; % set number of angle sampling points, e.g., from 0 to 180 degrees, one point every 2 degrees
theta_list = linspace(0, pi, N_angles); % angle list (radians)
fwhm_values = nan(1, N_angles); % store FWHM values calculated for each angle

% Prepare OTF needed for frequency domain convolution (already calculated before)
sheet_thickness = distance_per_pixel; % use thickness of a single pixel

% --- 2. Loop to calculate FWHM for each angle ---
% Create a timer to display estimated remaining time, as this part has heavy computation
h_waitbar = waitbar(0, 'Calculating FWHM for each angle...');
tic; % 启动计时器

for k = 1:N_angles
    % 更新等待条
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
    
    % a) 定义当前角度的法向量并创建3D薄片
    % 法向量定义为与z轴夹角为theta，在x-z平面内
    normal_vec = [sin(theta), 0, cos(theta)]; 
    plane_dist = X_real * normal_vec(1) + Y_real * normal_vec(2) + Z_real * normal_vec(3);
    sheet3D_slanted = double(abs(plane_dist) < (sheet_thickness / 2));
    
    % b) 频域方法计算卷积（成像）
    FT_sheet_slanted = fftn(fftshift(sheet3D_slanted));
    I_freq_slanted   = ifftshift(ifftn( FT_sheet_slanted .* OTF ));
    
    % c) 沿着法线方向提取1D强度剖面
    % 定义采样线
    max_dist = max(abs(dist_um)) * 0.8; % 在80%的范围内采样以避免边缘效应
    num_points = 1000; % 高密度采样以获得平滑曲线
    line_coords_1d = linspace(-max_dist, max_dist, num_points);
    % 计算采样点在3D空间中的坐标
    sample_points_x = line_coords_1d * normal_vec(1);
    sample_points_y = line_coords_1d * normal_vec(2);
    sample_points_z = line_coords_1d * normal_vec(3);
    % 使用三次插值提取强度值
    profile_1D = interp3(X_real, Y_real, Z_real, real(I_freq_slanted), ...
                           sample_points_x, sample_points_y, sample_points_z, 'cubic', 0);
                       
    % d) 计算剖面的FWHM
    if max(profile_1D) > 1e-9 % 确保剖面有有效的信号
        profile_norm = profile_1D / max(profile_1D);
        [~, pk_idx] = max(profile_norm); % 找到峰值位置
        
        % 寻找左侧半高点
        half_idx_left = find(profile_norm(1:pk_idx) <= 0.5, 1, 'last');
        % 寻找右侧半高点
        half_idx_right = find(profile_norm(pk_idx:end) <= 0.5, 1, 'first') + pk_idx - 1;
        
        if ~isempty(half_idx_left) && ~isempty(half_idx_right) && half_idx_left > 1
            % 通过线性插值精确计算半高点对应的坐标
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
close(h_waitbar); % 关闭等待条
fprintf('FWHM calculation finished in %.2f seconds.\n', toc);

% --- 3. 可视化结果 ---

fprintf('Starting analysis of FWHM (Method 2: Direct 2D IPSF Analysis)...\n');
% --- 1. 准备2D数据 ---
IPSF_xz = squeeze(IPSF(centre_pt, :, :)); % 提取 x-z 切片
% x_coords 和 z_coords 已经在前面定义为 dist_um
% --- 2. 循环计算每个角度的50%强度半径 ---
N_angles_direct = 180;                   
theta_list_direct = linspace(0, pi, N_angles_direct);  
x_half = nan(1, N_angles_direct);
z_half = nan(1, N_angles_direct);
% meshgrid for interp2
[Z_grid, X_grid] = meshgrid(z_coords, x_coords);

for k = 1 : N_angles_direct
    theta = theta_list_direct(k);
    
    % a) 确定当前角度下不超出图像边界的最大采样半径 r_max
    % 注意：您的代码中sin/cos定义与标准极坐标略有不同，这里保留您的定义
    % theta=0 沿x轴, theta=pi/2 沿z轴
    if abs(cos(theta)) < 1e-8,  r_max_x = Inf; else, r_max_x = max(abs(x_coords)) / abs(cos(theta)); end
    if abs(sin(theta)) < 1e-8,  r_max_z = Inf; else, r_max_z = max(abs(z_coords)) / abs(sin(theta)); end
    r_max = min(r_max_x, r_max_z);  
    
    % b) 沿该角度创建采样点
    N_r_sample = 800;
    r_vec = linspace(0, r_max, N_r_sample);
    
    % c) 计算采样点坐标并用插值法获取强度值
    x_r = r_vec * cos(theta);
    z_r = r_vec * sin(theta);
    f_r = interp2(Z_grid, X_grid, IPSF_xz, z_r, x_r, 'spline', 0 ); 
  
    % d) 归一化 (假设峰值在r=0处)
    if f_r(1) > 1e-9, f_r = f_r / f_r(1); else, f_r = zeros(size(f_r)); end
    
    % e) 寻找强度首次低于0.5的点，并精确计算其半径
    idx_fall = find( f_r < 0.5, 1, 'first' );
    if ~isempty(idx_fall) && idx_fall > 1
        r1 = r_vec(idx_fall - 1); r2 = r_vec(idx_fall);
        f1 = f_r(idx_fall - 1);   f2 = f_r(idx_fall);
        r_half = r1 + (0.5 - f1) * (r2 - r1) / (f2 - f1);
    else
        r_half = NaN;
    end
    
    % f) 将半径转换为轮廓点坐标
    x_half(k) =  r_half * cos(theta);
    z_half(k) =  r_half * sin(theta);
end
% --- 3. 创建完整轮廓线 ---
x_direct_full = [ x_half,  -fliplr(x_half) ];
z_direct_full = [ z_half,  -fliplr(z_half) ];
fprintf('Method 2 (Direct Analysis) complete.\n');
% a) 绘制 FWHM vs. Angle 的笛卡尔坐标图（信息最直接）
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
%ylim([0, max(fwhm_values)*1.1]); % 自动调整y轴范围

% b) 创建并绘制 FWHM 极坐标轮廓图，并以IPSF切片为背景
% 计算轮廓线的坐标
r_half = fwhm_values / 2;
x_contour = r_half .* sin(theta_list);
z_contour = r_half .* cos(theta_list);
x_contour_full = [x_contour, -fliplr(x_contour)];
z_contour_full = [z_contour, -fliplr(z_contour)];

% --- 开始绘图 ---
figure;
hold on; 

% 1. 绘制背景: IPSF 的 X-Z 切片
imagesc(z_coords, x_coords, IPSF_xz');
axis xy; 
colormap('parula');
hcb = colorbar;
hcb.Label.String = 'Normalized IPSF Intensity';
set(gca,'Color','k'); 

% 2. 叠加绘制第一条轮廓线 (Sheet FWHM，红色)
plot(x_contour_full, z_contour_full, 'g--', 'LineWidth', 2);

% 3. 叠加绘制第二条轮廓线 (IPSF FWHM，绿色)
plot(x_direct_full, z_direct_full, 'r--', 'LineWidth', 2); % 使用绿色虚线以示区别

% 4. 添加图表装饰和标签
set(gca, 'FontName', 'Times New Roman','FontSize', 14);
title('Comparison of FWHM Contours on IPSF', 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('z-axis (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('x-axis (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
lgd = legend('show', 'Location', 'best');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% 5. 设定显示范围和纵横比
daspect([1 1 1]); 
axis_limit_z = max(abs([z_contour_full, z_direct_full])) * 1.5;
axis_limit_x = max(abs([x_contour_full, x_direct_full])) * 1.5;
xlim([-1.5*zz,1.5*zz]);
ylim([-1.5*zz,1.5*zz]);
hold off;
fprintf('Final comparison visualization has been generated.\n');


%% Line Object (HORIZONTAL) Figures
%  ========================================================================
fprintf('\n--- LSF分析 PART 1: 沿Y轴的固定直线 ---\n');

% --- a) 创建一个沿X轴的理想直线对象 ---
line3D_y_axis = zeros(size(IPSF));
line3D_y_axis(centre_pt,: , centre_pt) = 1; % 

% --- b) 可视化直线对象 ---
figure;
p = patch(isosurface(X_real, Y_real, Z_real, line3D_y_axis, 0.5));
isonormals(X_real, Y_real, Z_real, line3D_y_axis, p);
p.FaceColor = 'red'; p.EdgeColor = 'none';
daspect([1,1,1]); view(3); axis tight; camlight; lighting gouraud;
grid on;
xlabel('x (µm)'); ylabel('y (µm)'); zlabel('z (µm)');
title('3D Line Object (along X-axis)');

% --- c) 频域卷积计算成像结果 (LSF) ---
FT_line = fftn(fftshift(line3D_y_axis));
OTF = ifftshift(fftn(fftshift(IPSF))); % calculate the OTF
OTF    = ifftshift(OTF);
I_LSF_3D = ifftshift(ifftn(FT_line .* OTF));
I_LSF_3D = real(I_LSF_3D); % 取实部，忽略计算误差

% --- d) 可视化LSF在X-Z平面的2D分布 ---
% LSF是线的像，其模糊分布在垂直于线的平面上，即x-z平面
lsf_xz_slice = squeeze(I_LSF_3D(centre_pt, :, :)); % 取y=0处的x-z切片
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
plot(0, 0, 'w+', 'MarkerSize', 10, 'LineWidth', 2); % 标记中心
daspect([1 1 1]);

% --- e) 提取并分析LSF的1D剖面 ---
lsf_z_profile = lsf_xz_slice(centre_pt, :); % 沿z轴的剖面 (x=0)

% 计算 Z-Profile 的 FWHM
fwhm_z = NaN; % 默认为NaN
profile_norm_z = lsf_z_profile' / max(lsf_z_profile); % 转置为列向量
[~, pk_idx_z] = max(profile_norm_z);

% 寻找左侧半高点
idx_left_z = find(profile_norm_z(1:pk_idx_z) <= 0.5, 1, 'last');
% 寻找右侧半高点
idx_right_z = find(profile_norm_z(pk_idx_z:end) <= 0.5, 1, 'first') + pk_idx_z - 1;

if ~isempty(idx_left_z) && ~isempty(idx_right_z) && idx_left_z > 1 && idx_right_z < length(profile_norm_z)
    % 线性插值精确计算左坐标
    y1L=profile_norm_z(idx_left_z);   x1L=z_coords(idx_left_z);
    y2L=profile_norm_z(idx_left_z+1); x2L=z_coords(idx_left_z+1);
    coord_L_z = x1L + (0.5 - y1L) * (x2L - x1L) / (y2L - y1L);
    % 线性插值精确计算右坐标
    y1R=profile_norm_z(idx_right_z-1); x1R=z_coords(idx_right_z-1);
    y2R=profile_norm_z(idx_right_z);   x2R=z_coords(idx_right_z);
    coord_R_z = x1R + (0.5 - y1R) * (x2R - x1R) / (y2R - y1R);
    
    fwhm_z = coord_R_z - coord_L_z;
end

% -- END: 更正后的FWHM计算方法 --


fprintf('FWHM of LSF along z-axis: %.4f µm\n', fwhm_z);

figure;
plot(z_coords, lsf_z_profile, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('Z-Profile (FWHM=%.2f µm)', fwhm_z));
grid on;
xlabel('Distance from center (µm)');
ylabel('Normalized Intensity');
title('1D Profiles of the Line Spread Function (LSF)');
legend('show', 'Location', 'best');
xlim([-2, 2]);

% --- f)【新增】水平线频域可视化 (修正版) ---
figure;
% 准备MTF背景 (u_x-u_z 平面)
% 注意：为了得到u_x-u_z平面，我们需要squeeze Y轴(第一个维度)
% 但为了和您其他代码一致，我们继续用squeeze X轴，并修正ylabel
otf_slice_2d = squeeze(MTF(centre_pt, :, :)); % 提取 u_y-u_z 平面
imagesc(spatial_freq_coords, spatial_freq_coords, otf_slice_2d);
axis xy image; colormap('parula');
caxis([0, 0.1]); 
hcb = colorbar; hcb.Label.String = 'MTF Value';
hold on;

% 对于沿X轴的线，其频谱在ux-uz平面上的交集是uz轴，即一条水平线 ux=0
% line(x_coords, y_coords)
line(spatial_freq_coords, zeros(size(spatial_freq_coords)), 'Color', 'red', 'LineWidth', 2);

% 添加正确的标签
xlabel('Axial Spatial Frequency u_z (\mum^{-1})');
ylabel('Lateral Spatial Frequency u_y (\mum^{-1})'); % 对应squeeze(MTF(centre_pt,:,:))
title({'Frequency View for Horizontal Line (along x-axis)', '(Red Line: Line Spectrum, Color: MTF)'});
xlim([-2*fc, 2*fc]); ylim([-2*fc, 2*fc]);
set(gca, 'Layer', 'top');
hold off;

% === Line (horizontal): 1D MTF along u_z at u_x=0  ===
% 取你之前用于频域可视化的切片：otf_slice_2d = squeeze(MTF(centre_pt, :, :));
% 该矩阵的"行"为 u_z，"列"为 u_x（见你前面的 imagesc 标注）
otf_slice_2d = squeeze(MTF(centre_pt, :, :));
% 找到 u_x = 0 的列索引
[~, idx_ux0] = min(abs(spatial_freq_coords - 0));
mtf_uz = otf_slice_2d(:, idx_ux0);              % 1D MTF vs u_z
mtf_uz = mtf_uz / max(mtf_uz(:) + eps);
uz = spatial_freq_coords(:);

% 只用非负频率半轴做 cutoff / MTF50 / MTF10
i0 = ceil(numel(uz)/2);   % 0 频附近索引（你的坐标通常对称、fftshift 后中心在中点）
uz_pos = uz(i0:end);
mtf_pos = mtf_uz(i0:end);

% 实操：直接用清晰的顺序语句最不容易踩坑（强烈建议用这一版）：
idx = find(mtf_pos <= 0.50, 1, 'first');
if isempty(idx) || idx<2
    fc50_lsf = NaN;
else
    fc50_lsf = uz_pos(idx-1) + (0.50 - mtf_pos(idx-1)) * ...
              (uz_pos(idx) - uz_pos(idx-1)) / (mtf_pos(idx) - mtf_pos(idx-1) + eps);
end

idx = find(mtf_pos <= 0.10, 1, 'first');
if isempty(idx) || idx<2
    fc10_lsf = NaN;
else
    fc10_lsf = uz_pos(idx-1) + (0.10 - mtf_pos(idx-1)) * ...
              (uz_pos(idx) - uz_pos(idx-1)) / (mtf_pos(idx) - mtf_pos(idx-1) + eps);
end

idx = find(mtf_pos <= 1e-3, 1, 'first');
if isempty(idx) || idx<2
    fcTH_lsf = NaN;
else
    fcTH_lsf = uz_pos(idx-1) + (1e-3 - mtf_pos(idx-1)) * ...
              (uz_pos(idx) - uz_pos(idx-1)) / (mtf_pos(idx) - mtf_pos(idx-1) + eps);
end

% 输出（把 µm 写成 um 也行，避免非 ASCII 字符带来的编码问题）
fprintf('Line (horizontal) 1D MTF along u_z: fc=%.3f cyc/um | MTF50=%.3f | MTF10=%.3f\n', ...
        fcTH_lsf, fc50_lsf, fc10_lsf);

% 画线前先判断
if isfinite(fc50_lsf), xline(fc50_lsf,':'); end
if isfinite(fc10_lsf), xline(fc10_lsf,':'); end
if isfinite(fcTH_lsf), xline(fcTH_lsf,'--'); end


% 小图：1D MTF（含三条阈值线）
figure; plot(uz_pos, mtf_pos, 'LineWidth', 1.8); grid on; hold on;
yline(0.50,':'); yline(0.10,':');
xline(fc50_lsf,':'); xline(fc10_lsf,':'); xline(fcTH_lsf,'--');
xlabel('u_z (cyc/\mum)'); ylabel('MTF (u_x=0)');
% 标题里的数值也用 num2str 规避 NaN（显示 N/A）
title(sprintf('1D MTF (u_x=0): fc=%.3f cyc/um | MTF50=%.3f | MTF10=%.3f\n', ...
        fcTH_lsf, fc50_lsf, fc10_lsf));


%% Line Object (ROTATE) Figures

fprintf('\n--- LSF分析 PART 2: 分析单一旋转角度的直线 ---\n');

% --- a) 构造并可视化指定角度的直线 ---

phi_rot_deg = 0; % << 您可以在这里设置任意角度，例如 30°
phi_rot_rad = phi_rot_deg * pi/180;

% 定义旋转后的方向向量 (在x-z平面内，从z轴向x轴旋转)

line_direction = [cos(phi_rot_rad), 0, sin(phi_rot_rad)];

% 使用到直线距离的公式创建3D直线对象

p_dot_d = X_real * line_direction(1) + Y_real * line_direction(2) + Z_real * line_direction(3);
dist_sq_from_line = (X_real - p_dot_d*line_direction(1)).^2 + ...
(Y_real - p_dot_d*line_direction(2)).^2 + ...
(Z_real - p_dot_d*line_direction(3)).^2;
line_radius = distance_per_pixel;
line3D_rotated = double(dist_sq_from_line < line_radius^2);

% 可视化旋转后的直线

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

% 用箭头标出直线方向

L = max(dist_um) * 0.8;
quiver3(0, 0, 0, line_direction(1), line_direction(2), line_direction(3), L, ...
'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
hold off;

% --- b) 频域卷积计算成像结果 (LSF) ---

FT_line_rot = fftn(fftshift(line3D_rotated));
I_LSF_rot_3D = ifftshift(ifftn(FT_line_rot .* OTF));
I_LSF_rot_3D = real(I_LSF_rot_3D); % 取实部

% --- c) 提取并分析LSF的1D剖面 (垂直于线方向) ---
% 定义垂直于线的采样方向

perp_direction = [-sin(phi_rot_rad), 0, cos(phi_rot_rad)];
max_dist_profile = 5;
num_points_profile = 1001;
line_coords_1d = linspace(-max_dist_profile, max_dist_profile, num_points_profile);

% 计算采样点坐标

sample_points_x = line_coords_1d * perp_direction(1);
sample_points_y = line_coords_1d * perp_direction(2);
sample_points_z = line_coords_1d * perp_direction(3);

% 插值提取剖面

profile_1D = interp3(X_real, Y_real, Z_real, I_LSF_rot_3D, ...
sample_points_x, sample_points_y, sample_points_z, 'cubic', 0);
profile_1D_norm = profile_1D / max(profile_1D);
fwhm_val = NaN; % 初始化为NaN, 以防计算失败
[~, pk_idx] = max(profile_1D_norm); % 找到峰值位置
idx_left = find(profile_1D_norm(1:pk_idx) <= 0.5, 1, 'last');
idx_right = find(profile_1D_norm(pk_idx:end) <= 0.5, 1, 'first') + pk_idx - 1;
if ~isempty(idx_left) && ~isempty(idx_right) && idx_left > 1 && idx_right < num_points_profile

% 线性插值精确计算左坐标

y1L = profile_1D_norm(idx_left); x1L = line_coords_1d(idx_left);
y2L = profile_1D_norm(idx_left+1); x2L = line_coords_1d(idx_left+1);
coord_L = x1L + (0.5 - y1L) * (x2L - x1L) / (y2L - y1L);

% 线性插值精确计算右坐标

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

% 将FWHM值加入标题，并明确是垂直(Normal)方向
set(gca, 'FontName', 'Times New Roman','FontSize', 14);
title(sprintf('LSF Profile Normal to Line (Angle = %.1f°, FWHM = %.3f µm)', phi_rot_deg, fwhm_val), 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Distance Normal to Line Direction (µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Normalized Intensity', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

%xlim([-2, 2]);
ylim([0,1]);
hold off;


% --- f) 【新增】在旋转坐标系中可视化LSF ---
fprintf('Resampling LSF into rotated coordinate system...\n');

% 1. 定义新的坐标轴方向：u轴垂直于线，v轴平行于线
u_vec = perp_direction; % u_vec is the normal vector
v_vec = line_direction; % v_vec is the parallel vector

% 2. 创建新的2D采样网格的坐标
% u轴（横轴）是观察模糊的方向，范围可以小一些
u_coords = linspace(-5, 5, 201); 
% v轴（纵轴）是沿着线的方向，范围可以大一些
v_coords = linspace(-5, 5, 401); 
[U_grid, V_grid] = meshgrid(u_coords, v_coords);

% 3. 将2D网格点映射到3D空间坐标
% P_3d = u * u_vec + v * v_vec
sample_points_X = U_grid .* u_vec(1) + V_grid .* v_vec(1);
sample_points_Y = U_grid .* u_vec(2) + V_grid .* v_vec(2); %
sample_points_Z = U_grid .* u_vec(3) + V_grid .* v_vec(3);

% 4. 使用interp3进行三维重采样
resampled_LSF_2D = interp3(X_real, Y_real, Z_real, I_LSF_rot_3D, ...
                           sample_points_X, sample_points_Y, sample_points_Z, 'cubic', 0);

resampled_LSF_2D_norm = resampled_LSF_2D / max(resampled_LSF_2D(:));

% 5. 绘制结果
figure;
imagesc(u_coords, v_coords, resampled_LSF_2D_norm');
axis xy; % 保证y轴（v轴）方向向上
colormap('parula'); % 使用'hot' colormap 通常能更好地显示光斑
colorbar;
daspect([1 1 1]); % 保证坐标轴比例正确

% 添加清晰的标签
colorbar
set(gca, 'FontName', 'Times New Roman','FontSize', 14);
title(sprintf('LSF in Rotated Frame (Line Angle = %.1f°)', phi_rot_deg), 'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Distance Normal to Line (u-axis, µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Distance Parallel to Line (v-axis, µm)', 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlim([-5,5]);
ylabel('Distance Parallel to Line (v-axis, µm)');


% --- g)【新增】旋转线频域可视化 ---
figure;
% 准备MTF背景
otf_slice_2d = squeeze(MTF(centre_pt, :, :)); % 提取 u_x-u_z 平面
imagesc(spatial_freq_coords, spatial_freq_coords, otf_slice_2d');
axis xy image; colormap('parula');
caxis([0, 0.1]);
hcb = colorbar; hcb.Label.String = 'MTF Value';
hold on;

% Calculate spectral line direction (perpendicular to spatial domain line direction)
% Spatial domain line direction (x,z): [cos(phi), sin(phi)]
% Frequency domain spectral line direction (ux,uz): [-sin(phi), cos(phi)]
spec_line_direction = [-sin(phi_rot_rad), cos(phi_rot_rad)];

% Draw spectral line
t_range = max(abs(spatial_freq_coords));
t = linspace(-t_range, t_range, 1000);
u_x_spec_line = t * spec_line_direction(1); % 频域x分量
u_z_spec_line = t * spec_line_direction(2); % 频域z分量
plot(u_x_spec_line, u_z_spec_line, 'r--', 'LineWidth', 2, 'DisplayName', 'Line Spectrum');

% 添加标签和标题
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
fprintf('\n--- LSF分析 PART 3: 计算FWHM vs. 旋转角度 ---\n');

% --- a) 设置旋转参数 ---
N_angles = 91; % 
phi_list_deg = linspace(0, 180, N_angles);
phi_list_rad = phi_list_deg * pi/180;
fwhm_lsf_values = nan(1, N_angles);

% --- b) 循环计算每个角度下LSF的FWHM ---
h_waitbar = waitbar(0, 'Calculating LSF FWHM for each line angle...');
tic;
for k = 1:N_angles
    % 更新进度条
    waitbar(k/N_angles, h_waitbar, sprintf('Processing angle %.1f°', phi_list_deg(k)));
    
    phi_rot_rad = phi_list_rad(k);
    
    % 1. 创建旋转后的直线对象 (与Part 2逻辑完全一致)
    line_direction = [cos(phi_rot_rad), 0, sin(phi_rot_rad)];
    p_dot_d = X_real * line_direction(1) + Y_real * line_direction(2) + Z_real * line_direction(3);
    dist_sq_from_line = (X_real - p_dot_d*line_direction(1)).^2 + ...
                        (Y_real - p_dot_d*line_direction(2)).^2 + ...
                        (Z_real - p_dot_d*line_direction(3)).^2;
    line_radius = distance_per_pixel;
    line3D_rotated = double(dist_sq_from_line < line_radius^2);
    
    % 2. 频域卷积 (与Part 2逻辑完全一致)
    FT_line_rot = fftn(fftshift(line3D_rotated));
    I_LSF_rot_3D = ifftshift(ifftn(FT_line_rot .* OTF));
    
    % 3. 沿垂直于线的方向提取剖面 (与Part 2逻辑完全一致)
    perp_direction = [-sin(phi_rot_rad), 0, cos(phi_rot_rad)]; 
    max_dist_profile = 5;
    num_points_profile = 1001;
    line_coords_1d = linspace(-max_dist_profile, max_dist_profile, num_points_profile);
    sample_points_x = line_coords_1d * perp_direction(1);
    sample_points_y = line_coords_1d * perp_direction(2); % 始终为0
    sample_points_z = line_coords_1d * perp_direction(3);
    
    profile_1D = interp3(X_real, Y_real, Z_real, real(I_LSF_rot_3D), ...
                         sample_points_x, sample_points_y, sample_points_z, 'cubic', 0);
                     
    % 4. 计算FWHM (与Part 2逻辑完全一致)
    fwhm_val = NaN;
    if max(profile_1D) > 1e-9 % 确保有有效信号
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
    
    % 5. 存储当前角度的FWHM值
    fwhm_lsf_values(k) = fwhm_val;
end
close(h_waitbar);
fprintf('LSF FWHM calculation for all angles finished in %.2f seconds.\n', toc);

% 可视化结果: FWHM vs. Angle overlay figure 坐标图 ---
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

% --- d) 可视化结果: LSF FWHM 轮廓图叠加在 IPSF 上 ---
% FWHM是垂直于线方向的宽度，所以轮廓图上的点，其方向就是perp_direction
r_lsf = fwhm_lsf_values / 2;
x_contour = r_lsf .* (-sin(phi_list_rad));
z_contour = r_lsf .* (cos(phi_list_rad));

% 创建完整轮廓 (0-360度)
x_contour2_full = [x_contour, -fliplr(x_contour), -x_contour, fliplr(x_contour)];
z_contour2_full = [z_contour, fliplr(z_contour), -z_contour, -fliplr(z_contour)];

figure;
hold on;
% 绘制背景: IPSF 的 X-Z 切片
I_LSF_3D = I_LSF_3D./max(I_LSF_3D(:));
IPSF_xz = squeeze(IPSF(centre_pt, :, :));
%I_freq_slanted、I_LSF_3D
imagesc(z_coords, x_coords, IPSF_xz');
axis xy; colormap('parula'); hcb = colorbar;
clim([0,1])

% 绘制 LSF FWHM 轮廓线
plot(x_contour2_full, z_contour2_full, 'y--', 'LineWidth', 2);
plot(x_contour_full, z_contour_full, 'g--', 'LineWidth', 2);
plot(x_direct_full, z_direct_full, 'r--', 'LineWidth', 2);


% 添加图表装饰
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
hold off;

fprintf('\n>>> LSF分析全部完成 <<<\n');



%% === 组件2: Polar overlay: Point vs Sheet vs Line (FWHM vs angle from z-axis) ===
% 依赖变量（若存在则直接用）：
% Sheet: theta_list (从 z 轴量)           , fwhm_values
% Line : phi_list_rad (线方向，从 x 轴量) , fwhm_lsf_values  ——> 法向方向其实与"从 z 轴量"的角度相同，可直接用 phi_list_rad
% Point: theta_list_direct (从 x 轴量), r_half(半径=FWHM/2) ——> 可选

% -------- 环形平滑（可选） --------
w = 7;  % 5–11 之间选一个奇数
pad   = @(x)[x(end-w+1:end); x(:); x(1:w)];
unpad = @(x)x(w+1:end-w);
circ_smooth = @(x) (numel(x)>2) .* unpad(movmean(pad(x(:)), w, 'omitnan')) + (numel(x)<=2).*x(:); % 保持列向量

% ====================== Sheet（法向 FWHM_n） ======================
ang_sheet = theta_list(:);          % 已是"从 z 轴量"，范围 0..pi
rS = fwhm_values(:);
rS = circ_smooth(rS);               % 可注释掉以关闭平滑
% 用 NaN 断开闭合连线，避免穿过原点的直线
angS_full = [ang_sheet; ang_sheet+pi; NaN];
rS_full   = [rS;       rS;          NaN];

% ====================== Line（法向 FWHM_\perp） ======================
% 线方向向量为 [cos(phi),0,sin(phi)]，其法向为 [-sin(phi),0,cos(phi)]；
% 法向与 z 轴夹角 α 满足 cos(α)=cos(phi) → α=phi（在 0..pi 内恒等），
% 因此直接用 phi_list_rad 作为"从 z 轴量"的角度即可。
ang_line = phi_list_rad(:);         % 0..pi
rL = fwhm_lsf_values(:);
rL = circ_smooth(rL);
angL_full = [ang_line; ang_line+pi; NaN];
rL_full   = [rL;       rL;          NaN];

% ====================== Point（IPSF 半高等值 FWHM）【可选】 ======================

ang_point = (pi/2) - theta_list(:);     % 先把"从x轴量"转成"从z轴量"
ang_point = mod(ang_point, pi);         % 映射到 [0, pi)

rP = 2 * r_half_vec(:);                 % FWHM = 2*半径

% —— 1) 去掉 NaN，并按角度升序排序（避免乱序造成跨段连线）——
valid = isfinite(ang_point) & isfinite(rP);
ang_point = ang_point(valid); 
rP        = rP(valid);

[ang_point, idx] = sort(ang_point); 
rP = rP(idx);

% （可选）在排序后再做环形平滑
% rP = circ_smooth(rP);

% —— 2) 在"重复角"处插入 NaN，打断径向连线 —— 
dup = [false; abs(diff(ang_point)) < 1e-9];  % 相邻角度几乎相同
ang_plot = ang_point;  r_plot = rP;
ang_plot(dup) = NaN;   r_plot(dup) = NaN;

% —— 3) 扩展到 0..2π，并在半周之间再断开一次 —— 
angP_full = [ang_plot; NaN; ang_plot + pi; NaN];
rP_full   = [r_plot;  NaN; r_plot;        NaN];


% ====================== 作图 ======================
fig = figure('Name','Polar overlay: resolution anisotropy','Color','w');
pax = polaraxes(fig);  hold(pax,'on');

% 统一极轴语义：0°=axial(z)，90°=lateral(x)
pax.ThetaZeroLocation = 'top';
pax.ThetaDir          = 'counterclockwise';
pax.ThetaTick         = [0 90 180 270];
pax.ThetaTickLabel    = {'0° axial (z)','90° lateral (x)','180°','270°'};

% 逐条绘制（存在才画）
polarplot(pax, angP_full, rP_full, '--', 'LineWidth', 2, 'DisplayName','Point  FWHM (IPSF)');
polarplot(pax, angS_full, rS_full, '--', 'LineWidth', 2, 'DisplayName','Sheet  FWHM_n');
polarplot(pax, angL_full, rL_full, '--',  'LineWidth', 2, 'DisplayName','Line   FWHM_\perp');

% 半径范围
r_candidates = [rS(:); rL(:)];
if ~isempty(rP), r_candidates = [r_candidates; rP(:)]; end
rmax = max(r_candidates, [], 'omitnan');
if isfinite(rmax) && rmax>0, rlim(pax, [0, 1.1*rmax]); end

title(pax,'FWHM vs orientation (angle from z-axis)');
legend(pax,'Location','southoutside');
grid(pax,'on');





%% === 组件1: MPA Dashboard: Wide-field Microscopy Performance Atlas =====
% This section generates a standardized performance dashboard 
% integrating spatial, frequency, and energy domain metrics

fprintf('\n========== Generating MPA Dashboard ==========\n');

% Create figure with specific size for publication quality
fig_dash = figure('Name', 'MPA Dashboard - Wide-field Microscopy', ...
                  'Position', [100, 50, 1400, 900], ...
                  'Color', 'white');

% Use tiledlayout for professional subplot arrangement
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Overall title for the dashboard
title(t, sprintf('Wide-field Microscopy Performance Dashboard (NA=%.1f, λ=%dnm)', ...
                 NA, lambda*1e9), ...
      'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Arial');

% === Subplot 1: PSF Lateral Profile (X) ===
nexttile;
% Extract lateral profile through center
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
    FWHM_lat = xR - xL;
    
    % Mark FWHM boundaries
    plot([xL, xL], [0, 0.5], 'r-', 'LineWidth', 1.5);
    plot([xR, xR], [0, 0.5], 'r-', 'LineWidth', 1.5);
    plot([xL, xR], [0.5, 0.5], 'r-', 'LineWidth', 2);
    
    % Add FWHM annotation
    text(0, 0.3, sprintf('FWHM = %.3f μm', FWHM_lat), ...
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

% === Subplot 2: Sheet Spread Function (Axial Profile) ===
nexttile;
% Use the axial response from sheet imaging
if exist('ifq', 'var') && exist('zp', 'var')
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
    if exist('FWHM_axial_um', 'var')
        text(0, 0.3, sprintf('FWHM_z = %.3f μm', FWHM_axial_um), ...
             'HorizontalAlignment', 'center', ...
             'FontSize', 12, 'FontWeight', 'bold', ...
             'Color', 'r', 'BackgroundColor', 'white', ...
             'EdgeColor', 'r', 'LineWidth', 1);
    end
end

% Formatting
grid on; box on;
xlim([-10, 10]);
ylim([0, 1.05]);
xlabel('Axial Position z (μm)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Normalized Intensity', 'FontSize', 11, 'FontWeight', 'bold');
title('Axial Sheet Response (Optical Sectioning)', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontName', 'Arial');

% === Subplot 3: Radial MTF with Key Metrics ===
nexttile;
% Use radial MTF if available
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
    text(0.8, 0.7, text_str, ...
         'FontSize', 11, 'FontWeight', 'bold', ...
         'BackgroundColor', 'white', ...
         'EdgeColor', 'k', 'LineWidth', 1);
else
    % Fallback: compute simple radial average
    MTF_z0 = squeeze(MTF(:,:,centre_pt));
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
    
    plot(fr_bins, mtf_radial/max(mtf_radial), 'b-', 'LineWidth', 2);
    hold on;
    yline(0.5, 'g--', 'LineWidth', 1, 'Alpha', 0.7);
    yline(0.1, 'r--', 'LineWidth', 1, 'Alpha', 0.7);
end

% Formatting
grid on; box on;
xlim([0, 3]);
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
p_outer.FaceColor = 'blue';  % Light blue
p_outer.EdgeColor = 'none';
p_outer.FaceAlpha = 0.11;  % Very transparent

hold on;

% Create inner isosurface at half-maximum (HMEF boundary)
p_inner = patch(isosurface(X_real, Y_real, Z_real, IPSF, half_max_level));
isonormals(X_real, Y_real, Z_real, IPSF, p_inner);

% Set appearance for inner surface (HMEF region)
p_inner.FaceColor = 'red';  % Darker blue
p_inner.EdgeColor = 'none';
p_inner.FaceAlpha = 1;  % More opaque

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
text(4, 4, 0, text_3d, ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', 'white', ...
     'EdgeColor', [0.2, 0.5, 0.8], 'LineWidth', 1.5, ...
     'HorizontalAlignment', 'center');

% Add legend explaining the visualization
text(3, 3, 4, sprintf('Red: FWHM region\nBlue: Full PSF'), ...
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
% saveas(fig_dash, 'MPA_Dashboard_Widefield.png');
% saveas(fig_dash, 'MPA_Dashboard_Widefield.fig');
% print(fig_dash, 'MPA_Dashboard_Widefield', '-dpng', '-r300');  % High resolution

fprintf('MPA Dashboard generated successfully!\n');
fprintf('Key Performance Metrics:\n');
fprintf('  - Lateral FWHM: %.3f μm\n', FWHM_lateral_um);
fprintf('  - Axial FWHM: %.3f μm\n', FWHM_f);
fprintf('  - Cutoff frequency: %.2f cyc/μm\n', fc_thresh);
fprintf('  - HMEF: %.3f\n', HMEF);
fprintf('========================================\n');



%% === 新组件2: ARP: Anisotropic Resolution Profile ===
% This enhanced version provides clearer visualization and better aesthetics

% Check which microscopy type we're plotting
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
        rtick_step = 1;
    end
    
    pax.RLim = [0, 6];
    pax.RTick = 0:rtick_step:6;
    
    % Add units to radial labels
    pax.RTickLabel = arrayfun(@(x) sprintf('%.1f μm', x), pax.RTick, 'UniformOutput', false);
    pax.RAxisLocation = 0;
end


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

%% ==新增代码：保存ARP数据 ==
fprintf('\n>>> 正在保存Wide-field ARP数据...\n');

% 创建一个结构体来存储数据
WF_data = struct();
WF_data.ang_PSF = angP_full;
WF_data.r_PSF   = rP_full;
WF_data.ang_SSF = angS_full;
WF_data.r_SSF   = rS_full;
WF_data.ang_LSF = angL_full;
WF_data.r_LSF   = rL_full;

% 将结构体保存到.mat文件
save('ARP_data_Widefield.mat', 'WF_data');

fprintf('>>> 数据成功保存至 ARP_data_Widefield.mat\n');
%%  PART 4: 辅助函数 (HELPER FUNCTION)
% ========================================================================
function fwhm = calculate_fwhm_smooth(coords_um, profile)
    fwhm = NaN;
    if max(profile(:)) < 1e-12, return; end
    profile_norm = profile / max(profile(:));
    factor = 20;
    xp = linspace(coords_um(1), coords_um(end), numel(coords_um)*factor);
    yp = interp1(coords_um, profile_norm, xp, 'spline');
    [~, pk_idx] = max(yp);
    if pk_idx > 1 && pk_idx < length(yp)
        find_half_points = @(y_data, x_data, y_target) ...
            deal(interp1(y_data(1:pk_idx), x_data(1:pk_idx), y_target, 'linear', NaN), ...
                 interp1(y_data(pk_idx:end), x_data(pk_idx:end), y_target, 'linear', NaN));
        [xL, xR] = find_half_points(yp, xp, 0.5);
        if ~isnan(xL) && ~isnan(xR), fwhm = abs(xR - xL); end
    end
end