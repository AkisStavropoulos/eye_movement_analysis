%% plot joystick output for SF for each subject separately
% have the files in the current folder each time
%% Janna
data1 = ImportSMR('s01_sf_nfb_b01_21-Aug-2018_task001.smr');
data2 = ImportSMR('s01_sf_nfb_b02_21-Aug-2018_task002.smr');
data3 = ImportSMR('s01_sf_nfb_b03_21-Aug-2018_task003.smr');

figure;plot(data1(9).imp.adc,data1(10).imp.adc,'.k','MarkerSize',.5);xlabel('angular velocity (deg/s)');ylabel('velocity (cm/s)');title('coefficient = 0');suptitle('Janna')
figure;plot(data2(9).imp.adc,data2(10).imp.adc,'.k','MarkerSize',.5);xlabel('angular velocity (deg/s)');ylabel('velocity (cm/s)');title('coefficient = 0.975');suptitle('Janna')
figure;plot(data3(9).imp.adc,data3(10).imp.adc,'.k','MarkerSize',.5);xlabel('angular velocity (deg/s)');ylabel('velocity (cm/s)');title('coefficient = 0.99');suptitle('Janna')

%% Henry
data1 = ImportSMR('s02_sf_nfb_b01_22-Aug-2018_task001.smr');
data2 = ImportSMR('s02_sf_nfb_b02_22-Aug-2018_task002.smr');
data3 = ImportSMR('s02_sf_nfb_b03_22-Aug-2018_task003.smr');

figure;plot(data1(9).imp.adc,data1(10).imp.adc,'.b','MarkerSize',.1);xlabel('angular velocity (deg/s)');ylabel('velocity (cm/s)');title('coefficient = 0');suptitle('Henry')
figure;plot(data2(9).imp.adc,data2(10).imp.adc,'.b','MarkerSize',.1);xlabel('angular velocity (deg/s)');ylabel('velocity (cm/s)');title('coefficient = 0.975');suptitle('Henry')
figure;plot(data3(9).imp.adc,data3(10).imp.adc,'.b','MarkerSize',.1);xlabel('angular velocity (deg/s)');ylabel('velocity (cm/s)');title('coefficient = 0.99');suptitle('Henry')

%% Baptiste
data1 = ImportSMR('s03_sf_nfb_b01_24-Aug-2018_task001.smr');
data2 = ImportSMR('s03_sf_nfb_b02_24-Aug-2018_task002.smr');
data3 = ImportSMR('s03_sf_nfb_b03_24-Aug-2018_task003.smr');

figure;plot(data1(9).imp.adc,data1(10).imp.adc,'.m','MarkerSize',.5);xlabel('angular velocity (deg/s)');ylabel('velocity (cm/s)');title('coefficient = 0');suptitle('Baptiste')
figure;plot(data2(9).imp.adc,data2(10).imp.adc,'.m','MarkerSize',.5);xlabel('angular velocity (deg/s)');ylabel('velocity (cm/s)');title('coefficient = 0.975');suptitle('Baptiste')
figure;plot(data3(9).imp.adc,data3(10).imp.adc,'.m','MarkerSize',.5);xlabel('angular velocity (deg/s)');ylabel('velocity (cm/s)');title('coefficient = 0.99');suptitle('Baptiste')
