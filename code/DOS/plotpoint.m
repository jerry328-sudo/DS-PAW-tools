for ii=1:76
    if dosdata.projectdos(ii).atom == 'O'
        scatter3(dosdata.projectdos(ii).position(1),dosdata.projectdos(ii).position(2),dosdata.projectdos(ii).position(3),'filled','r')
        hold on
    end
    if dosdata.projectdos(ii).atom == 'Ni'
        scatter3(dosdata.projectdos(ii).position(1),dosdata.projectdos(ii).position(2),dosdata.projectdos(ii).position(3), [],[0.9290 0.6940 0.1250],'filled')
        hold on
    end
    if dosdata.projectdos(ii).atom == 'C'
        scatter3(dosdata.projectdos(ii).position(1),dosdata.projectdos(ii).position(2),dosdata.projectdos(ii).position(3),'filled','b')
        hold on
    end
    if dosdata.projectdos(ii).atom == 'H'
        scatter3(dosdata.projectdos(ii).position(1),dosdata.projectdos(ii).position(2),dosdata.projectdos(ii).position(3),'filled','c')
        hold on
    end
    if dosdata.projectdos(ii).atom == 'Sn'
        scatter3(dosdata.projectdos(ii).position(1),dosdata.projectdos(ii).position(2),dosdata.projectdos(ii).position(3),'filled','g')
        hold on
    end
end
colormap hsv
axis equal