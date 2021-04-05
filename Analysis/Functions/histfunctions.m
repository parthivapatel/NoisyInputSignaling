hista = [];
figure()
hista = [hista; graphs(s005_x4)];
hista = [hista; graphs(s005_x5)];
hista = [hista; graphs(s005_x6)];
figure(10)
histogram(hista,20, 'FaceColor', [1 0 0], 'FaceAlpha', .2,'Normalization','probability')

histb = [];
figure()
histb = [histb; graphs(s01_x1)];
histb = [histb; graphs(s01_x2)];
histb = [histb; graphs(s01_x3)];
figure(10)
hold on
histogram(histb,20,'FaceColor', [0 1 0], 'FaceAlpha', .2,'Normalization','probability')

histc = [];
figure()
histc = [histc; graphs(s05_x2)];
histc = [histc; graphs(s05_x3)];
% graphs(s05_x1)
figure(10)
hold on
histogram(histc,20,'FaceColor', [0 0 1], 'FaceAlpha', .2,'Normalization','probability')

histd = [];
figure()
histd = [histd;graphs(s1_x4)];
histd = [histd;graphs(s1_x5)];
histd = [histd;graphs(s1_x6)];
figure(10)
hold on
histogram(histd,20,'FaceColor', [0 0 0], 'FaceAlpha', .2,'Normalization','probability')

histe = [];
figure()
histe = [histe;graphs(tnf1xy16)];
histe = [histe;graphs(tnf1xy17)];
histe = [histe;graphs(tnf1xy18)];
figure(10)
hold on
histogram(histe,20,'FaceColor', [0 0 0], 'FaceAlpha', .2,'Normalization','probability')

histf = [];
figure()
histf = [histf;graphs(tnf5xy4)];
histf = [histf;graphs(tnf5xy5)];
histf = [histf;graphs(tnf5xy6)];
figure(10)
hold on
histogram(histf,20,'FaceColor', [0 0 0], 'FaceAlpha', .2,'Normalization','probability')
