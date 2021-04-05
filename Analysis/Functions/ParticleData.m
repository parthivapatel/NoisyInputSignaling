function [particleNData] = ParticleData(tracks, k)
close all

particleNData = tracks(tracks(:,end)==k, :);

figure()
plot(particleNData(:,18),movmean(particleNData(:,5),2))

max(particleNData(:,5))
