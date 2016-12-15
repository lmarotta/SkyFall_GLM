%plot sample clips 
load Training_Data_labels.mat
ind1 = find(labels.value==1); %%Type of activity
clip = 121; %select a clip to plot
figure
plot(labels.acce{ind1(clip)}(:,2:end),'LineWidth',2)
xlabel('Sample')
ylabel('Acceleration [g]')
legend('X','Y','Z')
% save figure
ax = gca; 
ax.FontSize = 16;
fig = gcf;
print(fig,'./Figs/TripClip','-djpeg','-r300')