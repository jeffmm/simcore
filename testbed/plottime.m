% Plot the resultant time, speedup, and efficiency for the various
% versions of the force calculation

%brute, cells, cellsadj, neighbor list all pairs, neighbor cells
close all;

times = zeros(3,3,5);
times(:,:,1) = [2.04 1.52 1.19;
                85.16 43.29 24.82;
                1011.28 522.94 328.32];
times(:,:,2) = [2.42 1.70 1.31;
                22.04 11.33 6.42;
                13.66 8.42 6.04];
times(:,:,3) = [2.49 1.84 1.35;
                93.64 48.16 28.00;
                27.90 14.16 8.30];
times(:,:,4) = [1.42 1.16 1.01;
                9.63 5.18 3.08;
                27.17 15.30 9.12];
times(:,:,5) = [1.46 1.17 1.04;
                10.10 5.66 3.39;
                5.34 3.07 1.98];
ncores = [1 2 4];
natoms = [108 2916 78732];
scale = [1 10 500];
legendtext = {'108', '2916', '78732'};
xlegtext = 'ncores';
ylegtext = {'Time(s)', 'Speedup Sp', 'Efficiecy Ep'};
titletext = {'brute','cells','cells adj', 'neighbor all pairs', 'neighbor cells'};
% Now plot things
for i = 1:size(times,3)
    fig{i} = figure('Position',[100 100 1024 1200]);
    subplot(3,1,1);
    hold on;
    for j = 1:3
        plot(ncores,times(j,:,i),'LineWidth',2);
    end
    hold off;
    xlabel(xlegtext);
    ylabel(ylegtext(1));
    legend(legendtext,'FontSize',18);
    set(gca,'FontSize',18,'linewidth',2);
    title(titletext(i),'FontSize',18);
    
    % Calculate speedup Sp = Tp/Ts (time parallel over time serial)
    speedup = bsxfun(@ldivide,times(:,:,i),times(:,1,i));
    subplot(3,1,2);
    hold on;
    for j = 1:3
        plot(ncores,speedup(j,:),'LineWidth',2);
    end
    hold off;
    xlabel(xlegtext);
    ylabel(ylegtext(2));
    legend(legendtext,'FontSize',18);
    set(gca,'FontSize',18,'linewidth',2);
    
    % efficiency = Sp/p
    efficiency = bsxfun(@rdivide,speedup,ncores);
    subplot(3,1,3);
    hold on;
    for j = 1:3
        plot(ncores,efficiency(j,:),'LineWidth',2);
    end
    hold off;
    xlabel(xlegtext);
    ylabel(ylegtext(3));
    legend(legendtext,'FontSize',18);
    set(gca,'FontSize',18,'linewidth',2);
end

fig1 = figure('Position',[100 100 1024 1200]);
% Plot the 4 core version of each against natoms
hold on;
for i = 1:size(times,3)
    actualtime = times(:,3,i)'.*scale
    semilogy(actualtime);
end
hold off;