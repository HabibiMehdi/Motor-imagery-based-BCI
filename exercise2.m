

data = load('BCICIV_calib_ds1a_100Hz');

cnt = 0.1 * double(data.cnt);



chanlocs = repmat(struct('labels', [], 'X', [], 'Y', [], 'Z', [], 'sph_theta', [], 'sph_phi', [], 'theta', [], 'radius', [], 'type', []), 1, nChannels);
rotation_matrix = [0 1; 1 0];



% Fill in the chanlocs structure with your data
for i = 1:nChannels
    chanlocs(i).labels = clab{i};
    chanlocs(i).X = xpos(i);
    chanlocs(i).Y = ypos(i);
    chanlocs(i).Z = 0; % since this is a 2D projection, Z can be set to 0
    chanlocs(i).type = 'EEG';
    
    new_position = rotation_matrix * [chanlocs(i).X; chanlocs(i).Y];
    chanlocs(i).X = new_position(1);
    chanlocs(i).Y = new_position(2);

    % Convert Cartesian X, Y, Z to spherical and polar coordinates for topoplot
    [chanlocs(i).sph_theta, chanlocs(i).sph_phi, ~] = cart2sph(chanlocs(i).X, chanlocs(i).Y, chanlocs(i).Z);
    chanlocs(i).theta = rad2deg(chanlocs(i).sph_theta);
    chanlocs(i).radius = 0.6 * sqrt(chanlocs(i).X^2 + chanlocs(i).Y^2); 
end
radius = 0.6 * radius;
theta = theta + 270;

figure; 
topoplot([], chanlocs, 'style', 'blank', 'electrodes', 'labelpoint', 'chaninfo', 0.9);
title('EEG Channel Topography');

applyLaplacianAndPlot(chanlocs, cnt, false);


applyLaplacianAndPlot(chanlocs, cnt, true);

%%

plot(cnt(:,25));
title('plot ch25 before apply CAR');
%% apply CAR
mu = mean(cnt , 2); 

for i=1:size(cnt , 2)
    xi = cnt(:,i);
    cnt(:,i) = xi - mu ;
end
figure
plot(cnt(:,25),'r');
title('plot ch25 after apply CAR');