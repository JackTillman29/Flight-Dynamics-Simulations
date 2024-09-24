%% Detection Centroid (1-D)
% Find the centroids of connected cells in a random binary matrix
% detectionCentroid1d.m takes as input a binary detection vector and
% computes various parameters of each region (center, width)
%% Setup
% Create a smoothed gaussian noise signal
rng(128497)

N = [1 100];
x = abs(0.5*(randn(N) + 1i*randn(N))).^2;

kernel = ones(1,20);
kernel = kernel ./ sum(kernel(:));
PWR = conv(x,kernel,'same');
% Normalize
PWR = PWR ./ (1.1*max(PWR));

% Define a fixed threshold detector
det_thresh = 0.6;
A = PWR > det_thresh;
centroids = detectionCentroid1d(A,PWR);

% Isolate the parts of the signal over threshold
detPWR = PWR;
detPWR(A==0) = nan;

%% Plot
% Plot the original signal with portions above a fixed threshold colored
% red, and show the centroids of each region.
figure('Position',[520 378 889 420],'Color','w');
plot(PWR,'color',0.5*ones(1,3));
hold on;
plot(detPWR,'r.-','linewidth',3)
ylims = ylim;
plot(N,det_thresh*ones(1,2),'r--')
hc = stem(centroids.centroid_ind,ylims(2)+0*centroids.centroid_ind,'ro','MarkerFaceColor','y','LineWidth',2,'MarkerSize',10);
legend(hc,'Detection Region Centroid','location','best')
ylim(ylims+[0 0.1])
xlabel('Sample #')