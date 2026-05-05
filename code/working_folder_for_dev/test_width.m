function test(varargin)

p = inputParser;
addOptional(p, 'image_file_string', '../data/heart3003.nd2');
addOptional(p, 'min_rel_peak_prominence', 0.4);
addOptional(p, 'n_profile_points', 1000);

parse(p, varargin{:});
p = p.Results;

% Code

% Some marker points for spline
x = [433, 646. 761, 876, 921,  1036]
y = [198, 349, 567, 854, 1011, 1318]
% 
% x = [920, 962, 999, 1025, 1073, 1116, 1159, 1199]
% y = [1007, 1086, 1173, 1247, 1331, 1426, 1505, 1579];

% Fit the spline
cs = spline(x, y);
xs = linspace(x(1), x(end), 10000);
ys = ppval(cs, xs);

% Open the image
d = bfopen(p.image_file_string)

% Pull out first image
ds1 = d{1,1}{1,1};

% Extract profile along the spline
[prof_x, prof_y, im_profile] = improfile(ds1, xs, ys);

for i = 1 : (numel(prof_x) - 1)
    ang(i) = atan2((prof_y(i+1) - prof_y(i)), (prof_x(i+1) - prof_x(i)));
end

for i = 1 : (numel(prof_x) -1)
    % im_rot = rotateAround(ds1, prof_x(i), prof_y(i), -1 * ang(i) * 180 / pi);

    sx = 1;
    sy = 1;
    [H, W] = size(ds1);
    phi = ang(i);
    x0 = prof_x(i);
    y0 = prof_y(i);

    T = [sx*cos(phi), -sx*sin(phi), 0
         sy*sin(phi),  sy*cos(phi), 0
        (W+1)/2-((sx*x0*cos(phi))+(sy*y0*sin(phi))), (H+1)/2+((sx*x0*sin(phi))-(sy*y0*cos(phi))), 1];

    tform = affine2d(T);

    im_rot = imwarp(ds1, tform, 'OutputView', imref2d([H, W]), 'Interp', 'cubic');

    cx = W/2;
    cy = H/2;


    if (mod(i,10) == 1)
    
    
        figure(2);
        clf
        subplot(1,2,1);
        hold on;
        imagesc(ds1);
        plot(x0, y0, 'ro')
        ylim([1 H]);
        xlim([1 W]);
        set(gca, 'YDir', 'reverse');

        subplot(1,2,2);
        hold on;
        imagesc(im_rot);
        plot(cx, cy, 'rd')
        ylim([1 H]);
        xlim([1 W]);
        set(gca, 'YDir', 'reverse');
    
        drawnow;
        
    end
end



x_profile = 1:numel(im_profile);

% Find peaks
[pks, locs] = findpeaks(-im_profile, ...
    'MinPeakProminence', p.min_rel_peak_prominence * range(im_profile))

% Find sarcomere profiles
x_sarc_profile = linspace(0, 1, p.n_profile_points);

% Pull off profiles and interpolate
no_of_sarcomeres = numel(locs) - 1;
for i = 1 : no_of_sarcomeres
    profile_indices = locs(i) : locs(i+1);
    x_temp = normalize(profile_indices, 'range');
    y_sarc_profile(i, :) = interp1(x_temp, im_profile(profile_indices), ...
        x_sarc_profile);
end

figure(1);
clf;

n_rows = 4;
n_cols = 2;

subplot(n_rows, n_cols, 1);
imagesc(ds1);
hold on;
plot(x, y, 'k-o')
plot(xs, ys,'m-')

subplot(n_rows, n_cols, n_cols+1);
hold on;
plot(ang, 'b-');

% subplot(n_rows, n_cols, n_cols+1);
% hold on;
% plot(im_profile, 'b-');
% plot(locs, -pks, 'go')
% 
% subplot(n_rows, n_cols, (2*n_cols)+1);
% hold on;
% plot(x_sarc_profile, y_sarc_profile);
% 
% subplot(n_rows, n_cols, (3*n_cols)+1);
% hold on;
% shadedErrorBar(x_sarc_profile, y_sarc_profile, {@mean, @std});



