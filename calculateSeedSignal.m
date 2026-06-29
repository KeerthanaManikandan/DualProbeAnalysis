function seedSig = calculateSeedSignal(green, clipMask, seed, seedRad,frames)
% Calculating the averaged seed signal over a few pixels
% choose seed location (gaussian average of a radius around chosen point)
% seedRad = 6; % 6 = 250um radius, 12 = 500um radius, 24 = 1000um radius
clear gaussianDisk green_Seed clipMask_Seed
if exist('seedRad','var') == 0; seedRad = 6; end
seedSig = zeros(size(frames,3),1);

gaussianDisk = fspecial('gaussian',2*seedRad+1,seedRad);
green_seed = green(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad,:);
clipMask_seed = clipMask(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad);
gaussianDisk = gaussianDisk .* double(clipMask_seed);
gaussianDisk = gaussianDisk ./ sum(gaussianDisk(:));

for x = 1:2*seedRad+1     % rows are observations, columns are variables
    for y = 1:2*seedRad+1
        seedSig = seedSig + gaussianDisk(x,y) * squeeze(frames(seed(2)+x-seedRad-1,seed(1)+y-seedRad-1,:));

    end
end
end