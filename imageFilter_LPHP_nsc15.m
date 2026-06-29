function [imOut, imageBigHP] = imageFilter_LPHP_nsc15(imageBig,LPk,HPk,mask,HP_SS)
%{
close all; clearvars

% %% testing inputs
% imageBig = double(imread(['C:\Users\nsc15\Documents\Data\Merele_SqM\',...
%     'LeftHemisphere\02_28_2017\run00\Results\Ttest_nsc15\Frame8-9\LP10_HP150_C1.5\',...
%     '150BiPulses_60microAmps - Blank_submap_unFilt.bmp']));
imageBig = double(imread('C:\Users\nsc15\Documents\Data\Merele_SqM\LeftHemisphere\02_28_2017\run00\green00.bmp'));

LPk = 10;
HPk = 200/4;
mask = imread('C:\Users\nsc15\Documents\Data\Merele_SqM\LeftHemisphere\02_28_2017\run00\clipMask00.bmp');
mask = mask(:,:,1)>0;
%}

if ~exist('mask','var')
    mask = logical(ones(size(imageBig)));
end

if isempty(mask)
    mask = logical(ones(size(imageBig)));
end

if ~exist('HP_SS','var')
    HP_SS = 1;
end

if isempty(HP_SS)
    HP_SS = 1;
end

%% low pass
% imageBigLP = imfilter(imageBig, fspecial('gaussian', LPk, round(LPk/2)),'replicate');

if LPk>0
    mask2 = imfilter(double(mask), fspecial('gaussian', LPk, round(LPk/2)));
    mask3 = mask2+double(mask2 == 0);
    A_masked = imfilter(imageBig.*mask, fspecial('gaussian', LPk, round(LPk/2))); 
    A_masked = A_masked./mask3;
    imageBigLP = A_masked.*double(mask2 ~= 0) + imageBig.*double(mask2 == 0);
    imageBigLP(mask==0) = imageBig(mask==0);
else
    imageBigLP = imageBig;
end


%% high pass
if HPk > 0
    HPk = round(HPk/(4*HP_SS));
    if rem(HPk,2)~=0
        HPk = HPk+1;
    end

    maskSmall = imresize(mask,0.25/HP_SS,'bilinear');
    imageHP = imageBigLP;
    imageHP(~mask) = NaN;
    imageHP = imresize(imageHP,0.25/HP_SS,'bilinear');
    imageSmall = imresize(imageBig,0.25/HP_SS,'bilinear');
    imageP = NaN(size(imageHP,1)+HPk, size(imageHP,2)+HPk);
    imageP((HPk/2+1):(end-HPk/2),(HPk/2+1):(end-HPk/2)) = imageHP;

    imageF = NaN(size(imageHP));

    % dispstat('','init');
    for x = 1:size(imageHP,1)
        xW = x:x+HPk;

        for y = 1:size(imageHP,2)
            if ~isnan(imageHP(x,y))
                yW = y:y+HPk;
                imWindow = imageP(xW,yW);

                imageF(x,y) = median(imWindow(~isnan(imWindow)));
            else
                imageF(x,y) = 0;
            end
        end

    %     if rem(x-1,10)==0
    %         dispstat([num2str(round(x/size(imageHP,1)*10000)/100),'% Done.']);
    %     end
    end
    % dispstat('100% Done.');

    imageF(imageF==0) = median(median(imageF(maskSmall)));
    imageBigHP = imresize(imageF,[size(imageBig)],'bilinear');
    imageBigHP(~mask) = 0;

    imOut = imageBigLP - imageBigHP;
    disp('');
else
    imOut = imageBigLP;
end

% end