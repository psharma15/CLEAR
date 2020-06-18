function visImg1(imgX,nVoxel,voxImCoord,xyzSlice,roomSize)

p = [2 1 3];
imgX = reshape(imgX,nVoxel);
imgX = permute(imgX,p);

x = reshape(voxImCoord(:,1),nVoxel(2),nVoxel(1),nVoxel(3));
y = reshape(voxImCoord(:,2),nVoxel(2),nVoxel(1),nVoxel(3));
z = reshape(voxImCoord(:,3),nVoxel(2),nVoxel(1),nVoxel(3));


figure
h = slice(x,y,z,imgX,...
    xyzSlice{1},xyzSlice{2},xyzSlice{3});
xlabel('x-axis','FontSize',14)
ylabel('y-axis','FontSize',14)
zlabel('z-axis','FontSize',14)
xlim([roomSize(1,1), roomSize(1,2)])
ylim([roomSize(2,1), roomSize(2,2)])
zlim([roomSize(3,1), roomSize(3,2)])

set(h, 'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp');
alpha('color')

end
