figure; hold on;
axis([-1.1 1.1 -1.1 1.1]);
plot(empirical(:,1),empirical(:,2),'.','MarkerSize',10);
plot(dirac(:,1),dirac(:,2),'*');
for i=1:10
    text(empirical(i,1),empirical(i,2),num2str(i-1),'Color','red','FontSize',14)
end

for i=1:2
    text(dirac(i,1),dirac(i,2),num2str(i-1),'Color','red','FontSize',14)
end

