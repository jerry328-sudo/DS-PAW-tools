clc,clear
jsond=loadjson('dos.json');
for ii= 1:length(jsond.AtomInfo.Atoms)
    dos(ii).atom = jsond.AtomInfo.Atoms(1, ii).Element;
    dos(ii).position = jsond.AtomInfo.Atoms(1, ii).Position;
    count_orbit = 1;
    for ii1 = 1:length(jsond.DosInfo.Spin1.ProjectDos)
        if jsond.DosInfo.Spin1.ProjectDos(1, ii1).AtomIndex == ii
            orbit1{count_orbit} = jsond.DosInfo.Spin1.ProjectDos(1, ii1).Contribution;
            count_orbit = count_orbit +1;
        end
    end
    clear count_orbit
    try
        dos(ii).orbit.x = orbit1{1, 1};
        dos(ii).orbit.py = orbit1{1, 2};
        dos(ii).orbit.pz = orbit1{1, 3};
        dos(ii).orbit.px = orbit1{1, 4};
        dos(ii).orbit.dxy = orbit1{1, 5};
        dos(ii).orbit.dyz = orbit1{1, 6};
        dos(ii).orbit.dz2 = orbit1{1, 7};
        dos(ii).orbit.dxz = orbit1{1, 8};
        dos(ii).orbit.dx2 = orbit1{1, 9};
        dos(ii).orbit.p = orbit1{1, 2}+orbit1{1, 3}+orbit1{1, 4};
        dos(ii).orbit.d = orbit1{1, 5}+orbit1{1, 6}+orbit1{1, 7}+orbit1{1, 8}+orbit1{1, 9};
    catch
        disp([num2str(ii),'Element: ', jsond.AtomInfo.Atoms{1, ii}.Element, "'s"," orbit"," < ","9"]);
    end
end
dosdata.dosenergy = jsond.DosInfo.DosEnergy;
dosdata.efermi = jsond.DosInfo.EFermi;
dosdata.pdos = jsond.DosInfo.Spin1.Dos;
dosdata.projectdos = dos;
dosdata.dosnum = jsond.DosInfo.NumberOfDos;
clear dos ii ii1 jsond orbit1

%% 计算D带中心

% 对所有的d轨道进行加和
dband = zeros(1,dosdata.dosnum);
for ii = 1:length(dosdata.projectdos)
    dband = dband + dosdata.projectdos(ii).orbit.d;
end
dbande = dband.*dosdata.dosenergy;

dband_int = 0;
dbande_int = 0;
% 对d带积分
for ii = 2:length(dband)
    dband_int = dband_int + (dband(ii)+dband(ii-1))/2*(dosdata.dosenergy(ii)-dosdata.dosenergy(ii-1)); 
end
% 对de带积分

for ii = 2:length(dband)
    dbande_int = dbande_int + (dbande(ii)+dbande(ii-1))/2*(dosdata.dosenergy(ii)-dosdata.dosenergy(ii-1)); 
end

% 求d带中心
D_center = dbande_int/dband_int;
disp(['D带中心为:',num2str(D_center)]);

plot(dosdata.dosenergy,smooth(dband),'LineWidth',0.5)

hold on
plot(dosdata.dosenergy,smooth(dbande),'color','r','LineWidth',0.5)
% plot([D_center,D_center],[-600,600],'color','magenta','LineWidth',2)
h_a = gca;
set(h_a, 'XAxisLocation', 'origin');
set(h_a, 'YAxisLocation', 'origin');
% h_a.LineWidth = 1.5
xlim([-10 3])

dosenergy=dosdata.dosenergy';
dosdband=smooth(dband);
% dosdband=dosdband';
dosdbande=smooth(dbande);
test1 =[dosenergy,dosdband,dosdbande];

% save 'DosData.mat' dosdata
% save 'dos.txt' -ascii test1