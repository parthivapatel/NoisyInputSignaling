function [tracks] = ImageAnalysis(varargin)
close all


%for 8-bit avi

%bf = VideoReader(filenameBF);
%nmarker = VideoReader(strcat(filenameN,'_avi.avi')); %use avi file for segmentation
%reporter = VideoReader(filenameR);
%receptor = VideoReader(filenameG);
%bfsize = size(read(bf), 4);

%nmarker = dir(fullfile(strcat('.\',filenameN),'*.tif'));
%nmarker = {nmarker.name}';


% Base_Folder = 'E:\17_725_3t3_.05TNF_.2noise\Split';
% cd(Base_Folder)
% 
% BF_Folder = 'E:\17_725_3t3_.05TNF_.2noise\Split\C1-BF';
% YFP_Folder = 'E:\17_725_3t3_.05TNF_.2noise\Split\C3-RFP';
% Dapi_Folder = 'E:\17_725_3t3_.05TNF_.2noise\Split\C2-GFP';
% 
% cd(Dapi_Folder);
% folderlist = ls();

numvarargs = length(varargin);
if numvarargs > 6
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 5 optional inputs');
end

% set defaults for optional inputs (mincentroids, framerangemin,
% framerangemax, outliermin, outliermax, rng seed, receptor std threshhold)
optargs = {20 1 18 .2 2.5 1 .5};
optargs(1:numvarargs) = varargin;
% Place optional args in memorable variable names
[mincentroids, framerangemin, framerangemax, outliermin, outliermax, seed, receptorthresh] = optargs{:};
OutlierTolerancePercentage = [outliermin, outliermax];

results = cell(1, framerangemax-framerangemin+1);
finaltracks = [];
%centroiddata = cell(1, 70);

% Iterate over all frames
% bfsize
count = 1;
for i = 18
    
%     cd(strcat(Dapi_Folder,'\', folderlist(i,:)));
%     Dapi_Blank = imread(strcat('blank','t',num2str(1),'xy', num2str(i-2,'%02d'),'c',num2str(3),'.tif'));
%     cd(strcat(BF_Folder,'\', folderlist(i,:)));
%     BF_Blank = imread(strcat('blank','t',num2str(1),'xy', num2str(i-2,'%02d'),'c',num2str(1),'.tif'));
%     cd(strcat(YFP_Folder,'\', folderlist(i,:)));
%     YFP_Blank = imread(strcat('blank','t',num2str(1),'xy', num2str(i-2,'%02d'),'c',num2str(2),'.tif'));
%     
%     Dapi_Blank = im2double(Dapi_Blank)/max(max(im2double(Dapi_Blank)));
%     BF_Blank = im2double(BF_Blank)/max(max(im2double(BF_Blank)));
%     YFP_Blank = im2double(YFP_Blank)/max(max(im2double(YFP_Blank)));
        
    for k = framerangemin:framerangemax

        framenum = num2str(k,'%02d');
        numcentroids = 0;

        frame_n = im2double(imread(strcat('20180829_O1_Imm_SGC_XY',  num2str(i,'%02d'),'_C2.tif'), k));%./Dapi_Blank;
        frame_bf = im2double(imread(strcat('20180829_O1_Imm_SGC_XY', num2str(i,'%02d'),'_C1.tif'), k));%./BF_Blank;
        frame_r = im2double(imread(strcat('20180829_O1_Imm_SGC_XY', num2str(i,'%02d'),'_C3.tif'),k));%./YFP_Blank;

%         frame_n = (frame_n - min(frame_n(:))) / (max(frame_n(:)) - min(frame_n(:)));
        frame_n = im2uint8(frame_n);

%         frame_bf = (frame_bf - min(frame_bf(:))) / (max(frame_bf(:)) - min(frame_bf(:)));
        frame_bf = im2uint8(frame_bf);

%         frame_r = (frame_r - min(frame_r(:))) / (max(frame_r(:)) - min(frame_r(:)));
        frame_r = im2uint8(frame_r);

        % Iterate until imageprocessing captures at least mincentroids
        % centroids
        %fail = 0;

        %while numcentroids < mincentroids

        %watershed process starts 
        % imcomplement(segmented_images{2}))
        testnucleus = frame_n;
%        testnucleus = imcomplement(testnucleus);

        se = strel('disk',2);
        testnucleus = imdilate(testnucleus, se);
        testnucleus = imerode(testnucleus, se);
        testnucleus = imdilate(testnucleus, se);
        p = testnucleus';
        b = p(:)';  
        Y = sort(b);
%         length(b)

        testnucleus(testnucleus<= Y(500000)) = 0;
        testnucleus = bpass(testnucleus, 2, 51);
        testnucleus = imdilate(testnucleus, se);
        testnucleus = imerode(testnucleus, se);
        testnucleus = imerode(testnucleus, se);
%         testnucleus = imcomplement(testnucleus);
%         figure();imshow(testnucleus,[])


        %    imshow(tophatFiltered, [])
        %arbitrary
        %median(testnucleus)+.5*std2(testnucleus)
    %    testnucleus(testnucleus<= 10000) = 0;  

    %        testnucleus(testnucleus>= mean(max(testnucleus))-4.4*std2(testnucleus)) = 65535;
        %figure;imshow(testnucleus)
        testnucleus = imcomplement(testnucleus);
%         figure;imshow(testnucleus)

    %     if fail == 1
    %         k
    %         %color clustering on nucleomarker to isolate and color classify
    %         %nucleus center and edge with two different colors (colors 4+ will
    %         %usually start splitting up the background)
    % 
    %         cform = makecform('srgb2lab');
    %         colormark_n = applycform(frame_n,cform);
    % 
    %         n_colored = double(colormark_n(:,:,2:3));
    %         nrows = size(n_colored,1);
    %         ncols = size(n_colored,2);
    %         n_colored = reshape(n_colored,nrows*ncols,2);
    % 
    %         nColors = 4;
    %         % repeat the clustering 10 times to avoid local minima
    %         rng(seed);
    %         [cluster_idx, ~] = kmeans(n_colored,nColors,'distance','sqEuclidean', ...
    %                                               'Replicates',50);
    % 
    %         pixel_labels = reshape(cluster_idx,nrows,ncols);
    % 
    %         % optional figure to show clustering
    %         %figure; imshow(pixel_labels,[]), title('image labeled by cluster index');
    %         segmented_images = cell(1,3);
    %         rgb_label = repmat(pixel_labels,[1 1 3]);
    % 
    %         for j = 1:nColors
    %             color = frame_n;
    %             color(rgb_label ~= j) = 0;
    %             segmented_images{j} = color;
    %         end
    % 
    % 
    % 
    %         colornucleus = segmented_images{2} + segmented_images{3};
    %         I = rgb2gray(imcomplement(colornucleus));
    %     end

    %     % 5 radius disk for image processing
    %     se = strel('disk', 5);
    % 
    %     Ie = imerode(I, se);
    %     Iobr = imreconstruct(Ie, I);
    % 
    %     Iobrd = imdilate(Iobr, se);
    %     Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    %     Iobrcbr = imcomplement(Iobrcbr);
    % 
    %     Iobrcbr = imcomplement(Iobrcbr);
    % 
        bw = im2bw(testnucleus);
        bw = imcomplement(bw);

%         D = bwdist(~bw);
%         D = -D;
%         D(~bw) = -Inf;
% 
%         % Overlay with labels
%         DL = watershed(D);
        DL = bwlabel(bw);
%         figure;imshow(frame_n, 'InitialMag', 'fit')
    % 
    %     % Make a truecolor all-green image.
    %         green = cat(3, zeros(size(frame_bf)), ones(size(frame_bf)), zeros(size(frame_bf)));
    %         hold on
    %         h = imshow(green);
    %         hold off
    %         set(h, 'AlphaData', DL)
    %         % Show watershed ridgelines
    %         bgm = DL == 0; 
    %         hold on
    %         imshow(bgm), title('Watershed ridge lines (bgm)')
    %         hold off
    %         Lrgb = label2rgb(DL);
    %         figure; imshow(Lrgb), title('Watershed transform of gradient magnitude (Lrgb)')
    %         hold on
    %         CircleObjectsInImage(DL, 'r');
    %         hold off
    %       %frame_n_overlay = imread(strcat(filenameN,'.tif'),k); %rgb2gray frame_n

        % build properties out of watershed labels, and if the number of
        % centroids is lower than the minimum, repeat. (could change this
        % to a tolerance using standard deviation....)
        objectProperties = regionprops(DL,frame_n, ...
                'Area', 'Centroid', 'PixelList', 'PixelValues', 'MaxIntensity');

        numcentroids = length(objectProperties);
    %     if(numcentroids < mincentroids)
    %         fail = fail + 1;
    %         if fail > 1
    %             seed = seed + 1;
    %             
    %         end
    %     end


        % Brightfield properties
        objectProperties_bf = regionprops(DL,frame_bf, ...
                    'Area', 'Centroid', 'PixelList', 'PixelValues', 'MaxIntensity','MajorAxisLength','MinorAxisLength','Eccentricity', 'EquivDiameter', 'EulerNumber','Orientation','BoundingBox');

        % Receptor properties
    %    frame_g_overlay = frame_g; % green channel
    %    objectProperties_receptor = regionprops(DL,frame_g_overlay, ...
    %                'Area', 'Centroid', 'PixelList', 'PixelValues', 'MaxIntensity', 'BoundingBox' );

        % Reporter properties
        frame_r_overlay = frame_r; % Red channel
        objectProperties_reporter = regionprops(DL,frame_r_overlay, ...
                    'Area', 'Centroid', 'PixelList', 'PixelValues', 'MaxIntensity' );

        % create a list of outliers based on given tolerances for nuclear area
        medianMaxIntensity = median( [ objectProperties.Area ] );
        outliers = [ objectProperties.Area ] > OutlierTolerancePercentage(2) * medianMaxIntensity | [ objectProperties.Area ] < OutlierTolerancePercentage(1) * medianMaxIntensity;

        % remove outliers from watershed and objectproperites        
        if( sum( outliers ) > 0 )
            DL = RemoveObjectsFromImage( DL, objectProperties( outliers ) );
            objectProperties( outliers ) = [];  % remove outliers
            objectProperties_reporter( outliers ) = [];
            objectProperties_bf( outliers ) = [];
    %        objectProperties_receptor( outliers ) = [];
        end

        % grayscale texture features
         texture_featurelist = [];
    % %    totalcellbw = zeros(size(frame_g_overlay));
         for ii = 1:length(objectProperties_bf)
    % %        % capture cell boundary
    % %        xpos = objectProperties_bf(ii).Centroid(1);
    % %        ypos = objectProperties_bf(ii).Centroid(2);
    % %        circleprofile = [];
    % %        for iii = 1:5:360
    % %            degreeprofile = improfile(frame_g_overlay,[xpos,30*cosd(iii)], [ypos, 30*sind(iii)], 30);
    % %            ix = find(degreeprofile > max(degreeprofile)-10,1, 'first'); % - receptorthresh*std(degreeprofile)
    % %            xbound = xpos + ix*1*cosd(iii);
    % %            ybound = ypos + ix*1*sind(iii);
    % %            circleprofile = [circleprofile;[xbound,ybound]];
    % %        end
    % %        BWcellbound = roipoly(frame_g_overlay,circleprofile(:,1),circleprofile(:,2));
    % %        totalcellbw = totalcellbw + BWcellbound;
    %         
            subImage = imcrop(frame_bf, objectProperties_bf(ii).BoundingBox);
            Hausdimslope = [hausDim(subImage)];
            otsuthresh = otsurec(subImage, 5);
            sftafeats = sfta(subImage, 2 );

            % different angles: 0,45,90,135
            glcms_offsets = [0 1; -1 1;-1 0;-1 -1];
            [glcms, ~] = graycomatrix(subImage,'Offset',glcms_offsets);
            glcms_statlist = [];

            % iterate through different gray-level co-occurrence matrix at
            % different angles
            for iii = 1:length(glcms_offsets)
                glcms_stats = GLCMFeatures(glcms(:,:,iii));
                glcms_stats = struct2cell(glcms_stats);
                glcms_stats = cell2mat(glcms_stats);
                glcms_statlist = [glcms_statlist, glcms_stats.'];

            end
            % glcmsstats(19 features others)x4,
            % 1hausdimslope,5otsuthresholds,9sftafeats
            texture_featurelist = [texture_featurelist; glcms_statlist, Hausdimslope, otsuthresh.', sftafeats];      
        end
    %    figure;imshow(totalcellbw)   

        % receptor texture features
    %    texture_featurelist_r = [];
    %    for ii = 1:length(objectProperties_receptor)
    %        subImage_r = imcrop(frame_bf, objectProperties_receptor(ii).BoundingBox);
    %        Hausdimslope_r = [hausDim(subImage_r)];
    %        otsuthresh_r = otsurec(subImage_r, 5);
    %        sftafeats_r = sfta(subImage_r, 2 );

            % different angles: 0,45,90,135
    %        glcms_offsets_r = [0 1; -1 1;-1 0;-1 -1];
    %        [glcms_r, ~] = graycomatrix(subImage_r,'Offset',glcms_offsets_r);
    %        glcms_statlist_r = [];

            % iterate through different gray-level co-occurrence matrix at
            % different angles
    %        for iii = 1:length(glcms_offsets)
    %            glcms_stats_r = GLCMFeatures(glcms_r(:,:,iii));
    %            glcms_stats_r = struct2cell(glcms_stats_r);
    %            glcms_stats_r = cell2mat(glcms_stats_r);
    %            glcms_statlist = [glcms_statlist_r, glcms_stats_r.'];

    %        end
            % glcmsstats(19 features others)x4,
            % 1hausdimslope,5otsuthresholds,9sftafeats

    %        texture_featurelist_r = [texture_featurelist_r; glcms_statlist_r, Hausdimslope_r, otsuthresh_r.',sftafeats_r];      
    %    end

        XY = [];
        for ii = 1:length(objectProperties)
            XY =[XY;i];
        end

        % list of centroids
        centroids = [];
        for ii = 1:length(objectProperties)
            centroids =[centroids; objectProperties(ii).Centroid(1), objectProperties(ii).Centroid(2)];
        end

        % Bf scalar ojbectproperties
        etcinfo = [];
        for ii = 1:length(objectProperties_bf)
            etcinfo =[etcinfo; objectProperties_bf(ii).MajorAxisLength,objectProperties_bf(ii).MinorAxisLength,objectProperties_bf(ii).Eccentricity,...
                objectProperties_bf(ii).EquivDiameter, objectProperties_bf(ii).EulerNumber, objectProperties_bf(ii).Area];
        end

        % add framelabel
        numberOfCentroids = size(centroids, 1);
        centroids = [ centroids, k * ones( numberOfCentroids, 1 ) ];

        % distribution details
        bf_dist = [];
        for ii = 1:length(objectProperties_reporter)
            bf_dist =[bf_dist; mean(double(objectProperties_bf(ii).PixelValues)), std(double(objectProperties_bf(ii).PixelValues)),...
                sum(double(objectProperties_bf(ii).PixelValues))];
        end

        reporter_dist = [];
        for ii = 1:length(objectProperties_reporter)
            reporter_dist =[reporter_dist; mean(double(objectProperties_reporter(ii).PixelValues)), std(double(objectProperties_reporter(ii).PixelValues)), sum(double(objectProperties_reporter(ii).PixelValues))];
        end

        n_dist = [];
        for ii = 1:length(objectProperties_reporter)
            n_dist =[n_dist; mean(double(objectProperties(ii).PixelValues)), std(double(objectProperties(ii).PixelValues)), sum(double(objectProperties(ii).PixelValues))];
        end

    %    receptor_dist = [];
    %    for ii = 1:length(objectProperties_receptor)
    %        receptor_dist =[receptor_dist; mean(double(objectProperties_receptor(ii).PixelValues)), std(double(objectProperties_receptor(ii).PixelValues)), sum(double(objectProperties_receptor(ii).PixelValues))];
    %    end

        % find angle of ellipse orientation for each nucleus
        orientations = [];
        for ii = 1:length(objectProperties)
            orientations =[orientations; objectProperties_bf(ii).Orientation];
        end

%         % Image sift for each centroid (assume each nucleus is same scale and orientation for analysis purposes)
%         VL_data = [];
%         for ii = 1:length(centroids)
%             %orientations(ii) for orientation using ellipse angle instead of 0
%             [~,d] = vl_sift(single(frame_bf),'frames',[centroids(ii,1:2).';1;0]);
%             VL_data = [VL_data; d.'];
%         end

        %receptor_dist,  texture_featurelist_r, texture_featurelist, double(VL_data)
        framedata = [centroids(:,1:2) reporter_dist n_dist bf_dist etcinfo texture_featurelist XY centroids(:,3)];

        % centroiddata{k} = centroids;
        results{k} = framedata;
        count = count +1;
    %{
         for ii = 125
                % initial guess for sigma based on area of bright spot
                maximumPixelValue = max( double(objectProperties(ii).PixelValues) );
                darkPixelValue = median( double(objectProperties(ii).PixelValues) );
                pixelCountAboveHalf = sum( double(objectProperties(ii).PixelValues) > .5 *  ( maximumPixelValue + darkPixelValue ) );
                sigmaInitialGuess = 0.8 * sqrt( pixelCountAboveHalf / 2 / pi / log(2) );

                initialGuesses = [ ...
                    objectProperties(ii).Centroid(1), ... % yCenter
                    objectProperties(ii).Centroid(2), ... % xCenter
                    max(double(objectProperties(ii).PixelValues)) - min(double(objectProperties(ii).PixelValues)), ... % amplitude
                    sigmaInitialGuess, ... % (objectProperties(ii).BoundingBox(3) - 6) / 4, ... % sigma
                    min(double(objectProperties(ii).PixelValues)) ];

                BestFitData{ii} = nlinfit( objectProperties(ii).PixelList, double(objectProperties(ii).PixelValues), @Gaussian2DFitFunction, initialGuesses );

                % plot data, initial guess, and fit for each peak
                figure(17)
                clf

                % generate a triangle mesh from the best fit solution found by 
                % nlinfit and plot it
                gd = delaunay( objectProperties(ii).PixelList(:,1), ...
                    objectProperties(ii).PixelList(:,2) );
                trimesh( gd, objectProperties(ii).PixelList(:,1), ...
                    objectProperties(ii).PixelList(:,2), ...
                    Gaussian2DFitFunction(BestFitData{ii}, ...
                    objectProperties(ii).PixelList ) )
                hold on

                % plot initial guesses -- commented out to make plots less
                % cluttered. put this back in to debug initial guesses
                plot3( objectProperties(ii).PixelList(:,1), ...
                   objectProperties(ii).PixelList(:,2), ...
                   Gaussian2DFitFunction(initialGuesses, ...
                   objectProperties(ii).PixelList ), 'rx' )

                % plot image data
                plot3( objectProperties(ii).PixelList(:,1), ...
                    objectProperties(ii).PixelList(:,2), ...
                    objectProperties(ii).PixelValues, 'gx', 'LineWidth', 3)
                title(['Image data vs. Best Fit for Object Number ' num2str(ii)]);
        %}
    end
allData = vertcat( results{:} );
%centroiddata = vertcat(centroiddata{:});



% Use track program to track centroids
param=struct(); %structure array
param.mem = 4; %can leave for 10 frames
param.dim = 2;
param.good = 15; %set high when code is working (max is number of frames)
param.quiet = 1;
maxdisplacement = 15;
tracks = track(allData,maxdisplacement,param);


% colors = ['r','m','c','y']; 

% n = 1;
% particleNCentroids = tracks(tracks(:,4) == n, : );

% Trajectory1=particleNCentroids(:,1:2) * (7.4)/40;
% TauValues1=particleNCentroids(:,3);
% 
% Msd = MeanSquaredDisplacement( Trajectory1, TauValues1 );
% loglog(TauValues1*timescale,Msd, colors(n))
% hold on
% 
% n = 3;
% particleNCentroids = tracks(tracks(:,4) == n, : );
% 
% Trajectory3=particleNCentroids(:,1:2) * (7.4)/40;
% TauValues3=particleNCentroids(:,3);
% 
% Msd = MeanSquaredDisplacement( Trajectory3, TauValues3 );
% loglog(TauValues3*timescale,Msd, colors(n))
% hold on    
% 
% particleNCentroids = tracks(tracks(:,4) == 1, : );
% particleNplus1Centroids = tracks(tracks(:,4) == 3, : );
% 
% Trajectory = abs(particleNCentroids(:,1:2)* (7.4)/40 - particleNplus1Centroids(:,1:2) * (7.4)/40)/1.41;
% TauValues= particleNCentroids(:,3);
% Msd = MeanSquaredDisplacement( Trajectory, TauValues );
% 
% 
% %    alpha = 
% %    gstar = (2*kb*293.15)/(3*pi*(.84*10.^-6)*i*(2*pi/timestep)*Msd)
% %    gprime = abs(gstar)*cos(pi*
% 
% loglog(TauValues*timescale,Msd,'b')
% hold off
% title('MSD vs. Tau')
% xlabel('Log(Tau) (sec)')
% ylabel('Log(MSD) (µm^2)')
% legend('1','3', 'diff', 'Location', 'Northwest' ) 
% 
% 
% figure(2)
% plot(Trajectory1(:,2), Trajectory1(:,1),'b')
% title('Trajectory for Bead 1')
% xlabel('X Position (µm)')
% ylabel('Y Position (µm)') 
% 
% figure(3)
% plot(Trajectory3(:,2), Trajectory3(:,1),'r')
% title('Trajectory for Bead 3')
% xlabel('X Position (µm)')
% ylabel('Y Position (µm)') 
% 
% figure(4)
% plot(Trajectory(:,2),Trajectory(:,1),'g')
% title('Trajectory for Beads 3 - 1')
% xlabel('X Position (µm)')
% ylabel('Y Position (µm)')      
end

end

function CircleObjectsInImage( LabelImage, BorderColor )
    boundaries = bwboundaries( LabelImage );	
    numberOfBoundaries = size( boundaries );
    for k = 1 : numberOfBoundaries
        thisBoundary = boundaries{k};
        plot(thisBoundary(:,2), thisBoundary(:,1), 'Color', BorderColor, 'LineWidth', 2);
    end
end

function OutputBinaryImage = RemoveObjectsFromImage( InputBinaryImage, ObjectProperties )
    OutputBinaryImage = InputBinaryImage;
    eliminatedPixels = vertcat( ObjectProperties.PixelList );
    allObjectIndexes = sub2ind( size( InputBinaryImage ), ...
        eliminatedPixels(:, 2), eliminatedPixels(:,1) );
    OutputBinaryImage( allObjectIndexes ) = 0;
end

function out = Gaussian2DFitFunction( Parameters, Coordinates )
    yCenter = Parameters(1);
    xCenter = Parameters(2);
    amplitude = Parameters(3);
    sigma = Parameters(4);
    offset = Parameters(5);
    
    out = amplitude * ...
        exp( -(( Coordinates(:, 1) - yCenter ).^2 + ( Coordinates(:, 2) - xCenter ).^2 ) ...
        ./ (2 * sigma .^ 2 )) + offset;    
end

